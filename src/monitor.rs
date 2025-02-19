// a single radio monitoring a set frequency range

use crate::averagepsd::AveragePsd;
use nbfm_monitor_ui::messages::Message;

use soapysdr::Direction::Rx;
use soapysdr::Device;
use std::time::Instant;
use num_complex::Complex;

pub struct Monitor<const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> {
    samp_rate: usize,
    center_freq: usize,
    sdr_channel: usize,
    average_psd: AveragePsd<FFT_SIZE, AVG_WINDOW_SIZE>,
    dev: Device,
    sender: crossbeam_channel::Sender<Message>,
    recvr: crossbeam_channel::Receiver<Message>,
}

impl <const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
    pub fn new(dev: Device, send: crossbeam_channel::Sender<Message>, recv: crossbeam_channel::Receiver<Message>, samp_rate: usize, center_freq: usize) -> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
        // right now, only single-channel SDRs work
        let sdr_channel = 0;

        dev.set_frequency(Rx, sdr_channel, center_freq as f64, ()).expect("Failed to set frequency");

        println!("Setting sample rate to {}", samp_rate);
        dev.set_sample_rate(Rx, sdr_channel, samp_rate as f64).expect("Failed to set sample rate");

        let average_psd = AveragePsd::<FFT_SIZE, AVG_WINDOW_SIZE>::new(samp_rate, center_freq);

        return Monitor::<FFT_SIZE, AVG_WINDOW_SIZE> {
            samp_rate: samp_rate,
            sdr_channel: sdr_channel,
            dev: dev,
            average_psd: average_psd,
            center_freq: center_freq,
            sender: send,
            recvr: recv,
        }
    }


    pub fn start(&mut self) {
        println!("Starting stream...");
        let mut stream = self.dev.rx_stream::<Complex<f32>>(&[self.sdr_channel]).expect("Failed to open RX stream");
        println!("Fetching MTU...");
        let mtu: usize = stream.mtu().expect("Failed to get MTU");
        // the buffer we read IQ data into from the SDR
        let mut buf = vec![Complex::new(0.0, 0.0); mtu];

        stream.activate(None).expect("failed to activate stream");

        // throw away first MTU of samples while SDR boots/settles
        // see: https://pysdr.org/content/rtlsdr.html#gain-setting
        stream.read(&mut [&mut buf[..mtu]], 1_000_000).expect("error reading stream");

        // TODO fetch message(s) from crossbeam_channel, look for STOP message
        let mut should_continue = true;
        let read_size = buf.len();
        while should_continue {
            let len = stream.read(&mut [&mut buf[..read_size]], 1_000_000).expect("error reading stream");
            println!("Read {} bytes...", len);

            let frame = self.do_fft(len, &buf).expect("Could not do fft");
            self.sender.send(frame).expect("Cannot send to channel from thread");

            let msg = self.recvr.recv().expect("Cannot recv from channel in thread");
            if msg == (Message::Stop{}) {
                should_continue = false;
            }
        }
        stream.deactivate(None).expect("failed to deactivate stream");
    }

    fn do_fft(&mut self, len: usize, buf: &Vec<Complex<f32>>) -> Result<Message, &'static str> {

        // TODO:cleanup consider moving this into the struct
        // ignore some data so buf is a multiple of FFT_SIZE
        // can also pad with zeroes
        let sample_start_index = len % FFT_SIZE;

        // how many FFTs do we need to consume all the samples
        let fft_size_multiple = (len - sample_start_index) / FFT_SIZE;
        print!("Computing {} FFTs...", fft_size_multiple);
        let start = Instant::now();

        for i in 0..fft_size_multiple {
            let start = sample_start_index + (i * FFT_SIZE);
            // "samples" is clobbered by the FFT (computed in-place)
            let mut samples: [Complex<f32>; FFT_SIZE] = buf[start..start+FFT_SIZE].try_into().expect("incorrect sample length being passed to PSD function");
            self.average_psd.update(&mut samples);

        }
        let duration = start.elapsed();
        println!(" in {:?}", duration);
        // TODO:feature set up the actual "find n highest peaks" functionality
        // 1) find all signals above squelch threshold
        // 2) identify "peaks"? (could do n highest easily, or could deconflict
        //      peaks and then report all of them)
        // 3) push into idk stdout or something

        // 1. find max(psd)
        // 2. "zero out" a 12.5 kHz channel (nbfm) around the peak
        // 3. repeat
        // TODO keep freqs/values as i32 or change?
        let msg = Message::Frame {
            freqs: self.average_psd.get_psd().into_iter().map(|x| x as i32).collect(),
            values: self.average_psd.get_freq_range().into_iter().map(|x| x as i32).collect(),
        };

        return Ok(msg);
    }

    pub fn get_samp_rate(&self) -> usize {
        return self.samp_rate;
    }

    pub fn get_center_freq(&self) -> usize {
        return self.center_freq;
    }
}
