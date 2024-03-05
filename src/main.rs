use std::io::Write;
use std::fs::File;
use std::cmp::min;
use std::process;
use num_complex::Complex;
use std::sync::Arc;
use soapysdr::Direction::Rx;
use rustfft::{FftPlanner, Fft};
use rustfft::num_traits::Pow;
use simple_moving_average::{SMA, NoSumSMA};

// window of the running average
const WINDOW_SIZE: usize = 10;
// size of the fft
const FFT_SIZE: usize = 2048;
// sample rate of our SDR
const SAMP_RATE: usize = 3_200_000;
const FREQ: usize = 147_000_000;

fn main() {
    println!("Starting nbfm-listener...");

    let dev_filter = "driver=rtlsdr";
    let channel = 0;
    //let fname = "testfile";
    // how many samples to capture before closing the stream
    let mut num_samples = FFT_SIZE + 1;//i64::MAX;

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

    println!("{}", SAMP_RATE);
    dev.set_sample_rate(Rx, channel, SAMP_RATE as f64).expect("Failed to set sample rate");

    let mut stream = dev.rx_stream::<Complex<f32>>(&[channel]).unwrap();
    let mut buf = vec![Complex::new(0.0, 0.0); stream.mtu().unwrap()];

    //let mut outfile = BufWriter::new(File::create(fname).expect("error opening output file"));

    stream.activate(None).expect("failed to activate stream");

    while num_samples > 0 && !sb.caught() {
        let read_size = min(num_samples as usize, buf.len());
        let len = stream.read(&mut [&mut buf[..read_size]], 1_000_000).expect("error reading stream");
        //write_complex_file(&buf[..len], &mut outfile).unwrap();
        num_samples -= len
    }

    stream.deactivate(None).expect("failed to deactivate stream");

    // TODO move fft into the loop above for continuous operation

    // below this should move to a setup func
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(FFT_SIZE);

    let mut fft_averages = vec!(NoSumSMA::<f32, f32, WINDOW_SIZE>::new(); FFT_SIZE);

    // TODO can also use an integer multiple of FFT_SIZE
    // eg buf.len() % FFT_SIZE, then buf[remainder..] into the FFT
    let sample_start_index = buf.len() - FFT_SIZE;
    let mut fft_samp: [Complex<f32>; FFT_SIZE] = buf[sample_start_index..].try_into().unwrap();

    calc_psd(fft, &mut fft_samp, &mut fft_averages);
    // TODO figure out return from that func
}

fn calc_psd(fft: Arc<dyn Fft<f32>>, samples: &mut [Complex<f32>; FFT_SIZE], fft_averages: &mut Vec<NoSumSMA::<f32, f32, WINDOW_SIZE>>) {
    // https://pysdr.org/content/sampling.html#calculating-power-spectral-density

    // fft[0] -> DC
    // fft[len/2 + 1] -> nyquist f
    // 2 - len/2 are positive frequencies, each step is samp_rate/len(fft)
    // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/

    // TODO optimize using process_with_scratch()
    fft.process(samples);

    let fftbuf: Vec<f32> = samples.iter().map(|val| {
        val.norm().pow(2) / (FFT_SIZE * SAMP_RATE) as f32
    }).collect();

    let mut psd: Vec<f32> = fftbuf.iter().map(|val| {
        val.log10()
    }).collect();

    fftshift(&mut psd);
    // psd now contains our completed PSD calculation
    // next, set up our X axis (the frequency steps of the FFT)
    let mut freq_range: Vec<f32> = vec!(0.0; FFT_SIZE);
    for i in 0..FFT_SIZE {
        freq_range[i] = FREQ as f32 + (SAMP_RATE as f32/-2.0) + ((SAMP_RATE/FFT_SIZE) * i) as f32;
    }
    // TODO plot psd vs freq_range and... it should work?

    let mut buffer = File::create("fft-out").expect("cannot open output file");
    for elem in samples.iter() {
        write!(buffer, "{} ", elem).expect("error writing");
    }
    let mut freq_buffer = File::create("fft-freq").expect("cannot open output file");
    for elem in freq_range.iter() {
        write!(freq_buffer, "{} ", elem).expect("error writing");
    }

    // keep a running average of FFT values
    for (index, element) in psd.into_iter().enumerate() {
        fft_averages[index].add_sample(element);
    }
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

    // TODO holy crap write tests for this
}

// fn write_scalar_file<W: Write>(src_buf: &[f32], mut dest_file: W) -> io::Result<()> {
//     for sample in src_buf {
//         dest_file.write_f32::<LittleEndian>(*sample)?;
//     }
//     Ok(())
// }
//
// fn write_complex_file<W: Write>(src_buf: &[Complex<f32>], mut dest_file: W) -> io::Result<()> {
//     for sample in src_buf {
//         dest_file.write_f32::<BigEndian>(sample.re)?;
//         dest_file.write_f32::<LittleEndian>(sample.im)?;
//     }
//     Ok(())
// }
