// a single radio monitoring a set frequency range

use crate::averagepsd::AveragePsd;
use nbfm_monitor_ui::messages::Message;

use soapysdr::Direction::Rx;
use std::cmp::min;
use std::time::Instant;
use num_complex::Complex;

pub struct Monitor<const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> {
    samp_rate: usize,
    center_freq: usize,
    channel: usize,
    average_psd: AveragePsd<FFT_SIZE, AVG_WINDOW_SIZE>,
    dev: soapysdr::Device,
    num_samples: usize,
}

impl <const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
    pub fn new(samp_rate: usize, center_freq: usize, dev_filter: String) -> Monitor<FFT_SIZE, AVG_WINDOW_SIZE> {
        // right now, only single-channel SDRs work
        let channel = 0;

        // how many samples to capture before closing the stream
        //let mut num_samples = i64::MAX;
        let num_samples = samp_rate * 5; // sample for n seconds

        let devs = soapysdr::enumerate(&dev_filter[..]).expect("Error listing devices");
        let dev_args = match devs.len() {
            0 => {
                eprintln!("no matching SDR devices found");
                // TODO return an error
                panic!("can't find any SDRs");
            }
            1 => devs.into_iter().next().unwrap(),
            n => {
                eprintln!("{} devices found. Choose from:",n);
                for dev in devs {
                    eprintln!("\t'{}'", dev);
                }
                // TODO return an error
                panic!("multiple matching SDR devices");
            }
        };

        let dev = soapysdr::Device::new(dev_args).expect("Error opening device");

        dev.set_frequency(Rx, channel, center_freq as f64, ()).expect("Failed to set frequency");

        println!("Setting sample rate to {}", samp_rate);
        dev.set_sample_rate(Rx, channel, samp_rate as f64).expect("Failed to set sample rate");

        let average_psd = AveragePsd::<FFT_SIZE, AVG_WINDOW_SIZE>::new(samp_rate, center_freq);

        return Monitor::<FFT_SIZE, AVG_WINDOW_SIZE> {
            samp_rate: samp_rate,
            channel: channel,
            num_samples: num_samples,
            dev: dev,
            average_psd: average_psd,
            center_freq: center_freq,
        }
    }


    pub fn start(&mut self) {
        let sb = signalbool::SignalBool::new(
            &[signalbool::Signal::SIGINT], signalbool::Flag::Restart,
        ).unwrap();

        println!("Starting stream...");
        let mut stream = self.dev.rx_stream::<Complex<f32>>(&[self.channel]).expect("Failed to open RX stream");
        println!("Fetching MTU...");
        let mtu: usize = stream.mtu().expect("Failed to get MTU");
        // the buffer we read IQ data into from the SDR
        let mut buf = vec![Complex::new(0.0, 0.0); mtu];

        stream.activate(None).expect("failed to activate stream");

        // throw away first MTU of samples while SDR boots/settles
        // see: https://pysdr.org/content/rtlsdr.html#gain-setting
        stream.read(&mut [&mut buf[..mtu]], 1_000_000).expect("error reading stream");

        while self.num_samples > 0 && !sb.caught() {
            let read_size = min(self.num_samples as usize, buf.len());
            let len = stream.read(&mut [&mut buf[..read_size]], 1_000_000).expect("error reading stream");
            println!("Read {} bytes...", len);
            self.num_samples -= len;

            // TODO:cleanup consider moving this into the struct
            // ignore some data so buf is a multiple of FFT_SIZE
            // can also pad with zeroes
            let sample_start_index = len % FFT_SIZE;

            // how many FFTs do we need to consume all the samples
            // only compute every other FFT for performance reasons
            let fft_size_multiple = (len - sample_start_index) / FFT_SIZE;
            print!("Computing {} FFTs...", fft_size_multiple / 2);
            let start = Instant::now();

            for i in 0..fft_size_multiple {
                if i % 2 == 0 {
                    continue; // only use every other set of samples, for speed
                }
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
        }

        stream.deactivate(None).expect("failed to deactivate stream");

        // write out results
        // let mut outfile = File::create("fft-out").expect("cannot open output file");
        // let mut freq_outfile = File::create("fft-freq").expect("cannot open output file");

        // write the PSD values
        // for elem in self.average_psd.get_psd().iter() {
        //     writeln!(outfile, "{} ", elem).expect("error writing");
        // }
        //
        // // write the fft's frequency steps to plot against
        // for elem in self.average_psd.get_freq_range().iter() {
        //     writeln!(freq_outfile, "{} ", elem/1_000_000.).expect("error writing");
        // }

        // TODO keep freqs/values as i32 or change?
        let msg = Message::Frame {
            freqs: self.average_psd.get_psd().into_iter().map(|x| x as i32).collect(),
            values: self.average_psd.get_freq_range().into_iter().map(|x| x as i32).collect(),
        };
    }

    pub fn get_samp_rate(&self) -> usize {
        return self.samp_rate;
    }

    pub fn get_center_freq(&self) -> usize {
        return self.center_freq;
    }
}
