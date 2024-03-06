use std::io;
use std::f32::consts::PI;
use std::io::Write;
use std::fs::File;
use std::cmp::min;
use std::process;
use std::time::Instant;
use num_complex::Complex;
use std::sync::Arc;
use soapysdr::Direction::Rx;
use rustfft::{FftPlanner, Fft};
use rustfft::num_traits::Pow;
use simple_moving_average::{SMA, NoSumSMA};

// sample rate of our SDR
const SAMP_RATE: usize = 3_200_000; // 3.2 MHz
const FREQ: usize = 147_000_000; // where to tune the SDR
const FFT_SIZE: usize = 512;
// divide into 6.25kHZ channels
// NBFM is ~11 kHz signal so it'll get smeared but hopefully we'll be ok
// ... unfortunately, the channel width is also "FFTs per second" if using the full sample set
// (e.g. 10,000 FFT/s when using 10kHz channel width)

// window of the running average
const WINDOW_SIZE: usize = (SAMP_RATE / FFT_SIZE) * 5; // 10 seconds of averaging
                                                       // (only compute every *other* FFT right
                                                       // now)

fn main() {
    println!("Starting nbfm-listener...");

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

    dev.set_frequency(Rx, channel, FREQ as f64, ()).expect("Failed to set frequency");

    println!("Setting sample rate to {}", SAMP_RATE);
    dev.set_sample_rate(Rx, channel, SAMP_RATE as f64).expect("Failed to set sample rate");

    println!("Starting stream...");
    let mut stream = dev.rx_stream::<Complex<f32>>(&[channel]).expect("Failed to open RX stream");
    println!("Fetching MTU...");
    let mtu = stream.mtu().expect("Failed to get MTU");
    // the buffer we read IQ data into from the SDR
    let mut buf = vec![Complex::new(0.0, 0.0); mtu];

    //let mut outfile = BufWriter::new(File::create(fname).expect("error opening output file"));

    // this should move to a setup func or PSD struct or somethin
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(FFT_SIZE);

    let mut outfile = File::create("fft-out").expect("cannot open output file");
    let mut freq_outfile = File::create("fft-freq").expect("cannot open output file");

    let mut fft_averages = vec!(NoSumSMA::<f32, f32, WINDOW_SIZE>::new(); FFT_SIZE);

    stream.activate(None).expect("failed to activate stream");

    while num_samples > 0 && !sb.caught() {
        let read_size = min(num_samples as usize, buf.len());
        let len = stream.read(&mut [&mut buf[..read_size]], 1_000_000).expect("error reading stream");
        println!("Read {} bytes...", len);
        num_samples -= len;

        // ignore some data so buf is a multiple of FFT_SIZE
        let sample_start_index = len % FFT_SIZE;

        // how many FFTs do we need to consume all the samples
        let fft_size_multiple = (len - sample_start_index) / FFT_SIZE;
        print!("Computing {} FFTs...", fft_size_multiple / 2);
        let start = Instant::now();

        let mut psd: [f32; FFT_SIZE];
        for i in 0..fft_size_multiple {
            if i % 2 == 0 {
                continue; // only use every other set of samples, for speed
            }
            let start = sample_start_index + (i * FFT_SIZE);
            // NOTE "samples" is clobbered by the FFT (computed in-place)
            let mut samples: [Complex<f32>; FFT_SIZE] = buf[start..start+FFT_SIZE].try_into().expect("incorrect sample length being passed to PSD function");
            psd = calc_psd(&fft, &mut samples).expect("failed to calculate PSD");

            // keep a running average of FFT values
            for (index, element) in psd.into_iter().enumerate() {
                fft_averages[index].add_sample(element);
            }
        }
        let duration = start.elapsed();
        println!(" in {:?}", duration);
    }

    stream.deactivate(None).expect("failed to deactivate stream");

    // write the PSD values
    for elem in fft_averages.iter() {
        writeln!(outfile, "{} ", elem.get_average()).expect("error writing");
    }

    // set up our X axis (the frequency steps of the FFT)
    let mut freq_range: Vec<f32> = vec!(0.0; FFT_SIZE);
    for i in 0..FFT_SIZE {
        freq_range[i] = FREQ as f32 + (SAMP_RATE as f32/-2.0) + ((SAMP_RATE/FFT_SIZE) * i) as f32;
    }

    // write the fft's frequency steps to plot against
    for elem in freq_range.iter() {
        writeln!(freq_outfile, "{} ", elem).expect("error writing");
    }
}

// TODO don't use io::Error here, I guess, probably
fn calc_psd(fft: &Arc<dyn Fft<f32>>, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<[f32; FFT_SIZE]> {
    // https://pysdr.org/content/sampling.html#calculating-power-spectral-density

    // fft[0] -> DC
    // fft[len/2 + 1] -> nyquist f
    // 2 - len/2 are positive frequencies, each step is samp_rate/len(fft)
    // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/

    // both of these modify samples in-place
    hamming_window(samples).expect("failed to apply hamming window to samples");
    fft.process(samples);

    let mut psd: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
    for i in 0..FFT_SIZE {
        let sample: f32 = samples[i].norm().pow(2) / (FFT_SIZE * SAMP_RATE) as f32;
        psd[i] = 10.0 * sample.log10(); // convert to dB
    }

    fftshift(&mut psd.to_vec());
    // psd now contains our completed PSD calculation

    Ok(psd)
}

fn hamming_window(samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<()> {
    // https://math.stackexchange.com/questions/248849/hamming-window-understanding-formula
    // TODO this really oughta be computed and cached inside a struct for the fft
    let window: [f32; FFT_SIZE] = (0..FFT_SIZE).map(|n| 0.54 - 0.46 * ((2.0*PI*n as f32) / (FFT_SIZE as f32 - 1.0)).cos()).collect::<Vec<f32>>().try_into().expect("failed to compute hamming window");

    for i in 0..FFT_SIZE {
        samples[i] = samples[i].scale(window[i]);
    }

    Ok(())
}

// https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html
fn fftshift(fftbuf: &mut Vec<f32>) {
    let buflen: usize = fftbuf.len();
    assert!(buflen%2 == 0, "fftshift can only handle ffts with even length");
    // todo handle odd-length ffts, I guess, if I have to

    // take the back half of the array, reverse it, and put it on the front
    let neg_freqs = fftbuf.split_off(buflen/2);

    for elem in neg_freqs.into_iter().rev() {
        fftbuf.insert(0, elem);
    }
}

#[cfg(test)]
mod tests {
    use crate::fftshift;
    #[test]
    fn runs_at_all() {
        assert_eq!(1,1);
    }

    #[test]
    fn fftshift_works() {
        let mut startbuf: Vec<f32> = vec![0.0, 1.0, 2.0, 3.0, -4.0, -3.0, -2.0, -1.0];
        fftshift(&mut startbuf);
        let endbuf: Vec<f32> = vec![-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        assert_eq!(startbuf, endbuf);
    }
}
