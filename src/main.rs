mod averagepsd;

use std::io::Write;
use std::fs::File;
use std::cmp::min;
use std::process;
use std::time::Instant;
use num_complex::Complex;
use soapysdr::Direction::Rx;

//use piston_window::{EventLoop, PistonWindow, WindowSettings};
use plotters::prelude::*;
//use plotters_piston::{draw_piston_window};

use crate::averagepsd::AveragePsd;

// sample rate of our SDR
const SAMP_RATE: usize = 3_200_000; // 3.2 MHz
const CENTER_FREQ: usize = 147_000_000; // where to tune the SDR
const FFT_SIZE: usize = 512;
// divide into 6.25kHZ channels
// NBFM is ~11 kHz signal so it'll get smeared but hopefully we'll be ok
// ... unfortunately, the channel width is also "FFTs per second" if using the full sample set
// (e.g. 10,000 FFT/s when using 10kHz channel width)

// window of the running average
const AVG_WINDOW_SIZE: usize = (SAMP_RATE / FFT_SIZE) * 5; // 10 seconds of averaging
                                                       // (only compute every *other*
                                                       // FFT right now)

fn main() {
    let dev_filter = "driver=rtlsdr";
    let channel = 0;
    //let fname = "testfile";
    // how many samples to capture before closing the stream
    //let mut num_samples = i64::MAX;
    let mut num_samples = SAMP_RATE * 5; // sample for n seconds

    let devs = soapysdr::enumerate(&dev_filter[..]).expect("Error listing devices");
    let dev_args = match devs.len() {
        0 => {
            eprintln!("no matching devices found");
            process::exit(1);
        }
        1 => devs.into_iter().next().unwrap(),
        n => {
            eprintln!("{} devices found. Choose from:",n);
            for dev in devs {
                eprintln!("\t'{}'", dev);
            }
            process::exit(1);
        }
    };

    let sb = signalbool::SignalBool::new(
        &[signalbool::Signal::SIGINT], signalbool::Flag::Restart,
    ).unwrap();

    let dev = soapysdr::Device::new(dev_args).expect("Error opening device");

    dev.set_frequency(Rx, channel, CENTER_FREQ as f64, ()).expect("Failed to set frequency");

    println!("Setting sample rate to {}", SAMP_RATE);
    dev.set_sample_rate(Rx, channel, SAMP_RATE as f64).expect("Failed to set sample rate");

    println!("Starting stream...");
    let mut stream = dev.rx_stream::<Complex<f32>>(&[channel]).expect("Failed to open RX stream");
    println!("Fetching MTU...");
    let mtu = stream.mtu().expect("Failed to get MTU");
    // the buffer we read IQ data into from the SDR
    let mut buf = vec![Complex::new(0.0, 0.0); mtu];

    // this should move to a setup func or PSD struct or somethin
    let mut average_psd: AveragePsd<AVG_WINDOW_SIZE, FFT_SIZE> = AveragePsd::<AVG_WINDOW_SIZE, FFT_SIZE>::new(SAMP_RATE, CENTER_FREQ);

    stream.activate(None).expect("failed to activate stream");

    while num_samples > 0 && !sb.caught() {
        let read_size = min(num_samples as usize, buf.len());
        let len = stream.read(&mut [&mut buf[..read_size]], 1_000_000).expect("error reading stream");
        println!("Read {} bytes...", len);
        num_samples -= len;

        // TODO consider moving this into the struct
        // ignore some data so buf is a multiple of FFT_SIZE
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
            average_psd.update(&mut samples);

        }
        // TODO set up the actual "find n highest peaks" functionality
        // 1) find all signals above squelch threshold
        // 2) identify "peaks"? (could do n highest easily, or could deconflict
        //      peaks and then report all of them)
        // 3) push into idk stdout or something
        let duration = start.elapsed();
        println!(" in {:?}", duration);
    }

    stream.deactivate(None).expect("failed to deactivate stream");

    // display the output using plotters-rs

    let drawing_area = SVGBackend::new("chart_builder_on.svg", (600, 400)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();

    let frequency_range = average_psd.get_freq_range();
    let data = average_psd.get_psd();

    let mut cc = ChartBuilder::on(&drawing_area)
        .margin(10)
        .caption("PSD of Spectrum", ("sans-serif",30))
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(frequency_range[0]..frequency_range[FFT_SIZE-1],
                            -120.0 as f32..20.0 as f32) // dB
        .unwrap();

    cc.configure_mesh()
        .x_label_formatter(&|x| format!("{}", x/1e6)) // convert Hz -> MHz
        //.y_label_formatter(&|y| format!("{}%", (*y * 100.0) as u32))
        .x_labels(10)
        .y_labels(5)
        .x_desc("Frequency (MHz)")
        .y_desc("Power (dB)")
        .axis_desc_style(("sans-serif", 15))
        .draw()
        .unwrap();

    println!("{:?}", frequency_range.iter().zip(data.iter()));
    cc.draw_series(LineSeries::new(
            frequency_range.iter().zip(data.iter()).map(|(a, b)| (*a, *b)),
            BLACK)
        ).unwrap();

    // write out results
    let mut outfile = File::create("fft-out").expect("cannot open output file");
    let mut freq_outfile = File::create("fft-freq").expect("cannot open output file");

    // write the PSD values
    for elem in average_psd.get_psd().iter() {
        writeln!(outfile, "{} ", elem).expect("error writing");
    }

    // write the fft's frequency steps to plot against
    for elem in average_psd.get_freq_range().iter() {
        writeln!(freq_outfile, "{} ", elem).expect("error writing");
    }
}


#[cfg(test)]
mod tests {
    #[test]
    fn runs_at_all() {
        assert_eq!(1,1);
    }

}
