// a single radio monitoring a set frequency range

use crate::averagepsd::AveragePsd;
use nbfm_monitor_ui::messages::Message;

use soapysdr::Direction::Rx;
use soapysdr::Device;
use num_complex::Complex;
use crossbeam_channel::TryRecvError;
use log::{error, info, debug};

pub struct Monitor<const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> {
    sdr_channel: usize,
    average_psd: AveragePsd<FFT_SIZE, AVG_WINDOW_SIZE>,
    dev: Device,
    sender: crossbeam_channel::Sender<Message>,
    pub recvr: crossbeam_channel::Receiver<Message>,
}

impl <const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
    pub fn new(dev: Device, send: crossbeam_channel::Sender<Message>, recv: crossbeam_channel::Receiver<Message>, samp_rate: usize, center_freq: usize) -> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
        // right now, only single-channel SDRs work
        let sdr_channel = 0;

        dev.set_frequency(Rx, sdr_channel, center_freq as f64, ()).expect("Failed to set frequency");

        debug!("Sample rate set to {}", samp_rate);
        dev.set_sample_rate(Rx, sdr_channel, samp_rate as f64).expect("Failed to set sample rate");

        let average_psd = AveragePsd::<FFT_SIZE, AVG_WINDOW_SIZE>::new(samp_rate, center_freq);

        return Monitor::<FFT_SIZE, AVG_WINDOW_SIZE> {
            sdr_channel: sdr_channel,
            dev: dev,
            average_psd: average_psd,
            sender: send,
            recvr: recv,
        }
    }


    pub fn start(&mut self) {
        info!("Starting stream...");
        let mut stream = self.dev.rx_stream::<Complex<f32>>(&[self.sdr_channel]).expect("Failed to open RX stream");
        let mtu: usize = stream.mtu().expect("Failed to get MTU");
        debug!("MTU set to {mtu}");
        // the buffer we read IQ data into from the SDR
        let mut buf = vec![Complex::new(0.0, 0.0); mtu];

        stream.activate(None).expect("failed to activate stream");

        // throw away first MTU of samples while SDR boots/settles
        // see: https://pysdr.org/content/rtlsdr.html#gain-setting
        stream.read(&mut [&mut buf[..mtu]], 1_000_000).expect("error reading SDR stream");

        let read_size = buf.len();
        loop {
            // fetch message from crossbeam_channel, look for STOP
            match self.recvr.try_recv() {
                Ok(Message::Stop{}) => {
                    break;
                },
                Ok(_) => {}, // ignore other message types for now
                Err(TryRecvError::Empty) => {},
                Err(e) => panic!("Panic reading channel in monitor: {e:?}"),
            }
            let len = match stream.read(&mut [&mut buf[..read_size]], 1_000_000) {
                Ok(length) => length,
                Err(e) => match e.code {
                    soapysdr::ErrorCode::Overflow => {
                        error!("SDR read overflow");
                        continue;
                    },
                    _ => panic!("{}", e),
                }
            };
            //println!("Read {} bytes...", len);

            // TODO parallelize more so we get fewer/no overflows
            // https://docs.rs/rayon/latest/rayon/ maybe
            // parallelize each FFT within a buf probably best (vs starting the
            // next buf of data before this one has finished) - but async on
            // the send(frame) would probably help too
            let frame = self.do_fft(len, &buf).expect("Could not do fft");
            self.sender.send(frame).expect("Cannot send to channel from thread");
        }

        // cleanup
        stream.deactivate(None).expect("failed to deactivate stream");
    }

    fn do_fft(&mut self, len: usize, buf: &Vec<Complex<f32>>) -> Result<Message, &'static str> {
        // ignore some data so buf is a multiple of FFT_SIZE
        // can also pad with zeroes
        let sample_start_index = len % FFT_SIZE;

        // how many FFTs do we need to consume all the samples
        let fft_size_multiple = (len - sample_start_index) / FFT_SIZE;
        //print!("Computing {} FFTs...", fft_size_multiple);
        //let start = Instant::now();

        for i in 0..fft_size_multiple {
            let start = sample_start_index + (i * FFT_SIZE);
            // "samples" is clobbered by the FFT (computed in-place)
            let mut samples: [Complex<f32>; FFT_SIZE] = buf[start..start+FFT_SIZE].try_into().expect("incorrect sample length being passed to PSD function");
            self.average_psd.update(&mut samples);

        }
        //let duration = start.elapsed();
        //println!(" in {:?}", duration);
        // TODO:feature set up the actual "find n highest peaks" functionality
        // 1) find all signals above squelch threshold
        // 2) identify "peaks"? (could do n highest easily, or could deconflict
        //      peaks and then report all of them)
        // 3) push into idk stdout or something

        // 1. find max(psd)
        // 2. "zero out" a 12.5 kHz channel (nbfm) around the peak
        // 3. repeat
        let msg = Message::Frame {
            values: self.average_psd.get_psd().into_iter().map(|x| x as f32).collect(),
            freqs: self.average_psd.get_freq_range().into_iter().map(|x| x as i32).collect(),
        };

        return Ok(msg);
    }
}
