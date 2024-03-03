use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::cmp::min;
use std::process;
use std::i64;
use byteorder::{LittleEndian, WriteBytesExt};
use num_complex::Complex;
use soapysdr::Direction::Rx;
use rustfft::FftPlanner;
use simple_moving_average::{SMA, NoSumSMA};

fn main() {
    println!("Starting nbfm-listener...");

    let dev_filter = "driver=rtlsdr";
    let channel = 0;
    let fname = "testfile";
    let mut num_samples = i64::MAX;

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

    let freq = 144e6;
    dev.set_frequency(Rx, channel, freq, ()).expect("Failed to set frequency");

    const SAMP_RATE: usize = 3_200_000;
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
        num_samples -= len as i64
    }

    // TODO move fft into the loop above for continuous operation

    stream.deactivate(None).expect("failed to deactivate stream");
    println!("stream deactivated, starting FFT...");

    // fft[0] -> DC
    // fft[len/2 + 1] -> nyquist f
    // 2 - len/2 are positive frequencies, each step is samp_rate/len(fft)
    // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/
    let fft_size = buf.len();
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_size);

    // TODO optimize using process_with_scratch()
    fft.process(&mut buf);

    let fftbuf: Vec<f32> = buf.iter().map(|val| {val.norm()}).collect();

    let mut fftoutfile = BufWriter::new(File::create("fft-output").expect("error opening output file"));
    write_scalar_file(&fftbuf, &mut fftoutfile).expect("failed to write fft output");

    // keep a running average of FFT values
    const WINDOW_SIZE: usize = 10;
    let mut fft_averages = vec!(NoSumSMA::<f32, f32, WINDOW_SIZE>::new(); fft_size);

    for (index, element) in fftbuf.into_iter().enumerate() {
        fft_averages[index].add_sample(element);
    }
}

fn write_scalar_file<W: Write>(src_buf: &[f32], mut dest_file: W) -> io::Result<()> {
    for sample in src_buf {
        dest_file.write_f32::<LittleEndian>(*sample)?;
    }
    Ok(())
}

fn write_complex_file<W: Write>(src_buf: &[Complex<f32>], mut dest_file: W) -> io::Result<()> {
    for sample in src_buf {
        dest_file.write_f32::<LittleEndian>(sample.re)?;
        dest_file.write_f32::<LittleEndian>(sample.im)?;
    }
    Ok(())
}
