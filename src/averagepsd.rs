use std::f32::consts::PI;
use rustfft::{FftPlanner, Fft};
use simple_moving_average::{SMA, NoSumSMA};
use std::sync::Arc;
use num_complex::Complex;
use std::io;


// keeps a running average of the power spectral density of a frequency range
pub struct AveragePsd<const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> {
    samp_rate: usize,
    center_freq: usize,
    fft: Arc<dyn Fft<f32>>,
    psd_averages: Vec<NoSumSMA::<f32, f32, { AVG_WINDOW_SIZE }>>,
    // probably worth implementing this myself soonish so that the window
    // can be adjusted on the fly rather than passing a const
}

impl<const FFT_SIZE: usize, const AVG_WINDOW_SIZE: usize> AveragePsd<FFT_SIZE, AVG_WINDOW_SIZE> {
    //const PSD_SIZE: usize = (FFT_SIZE / 2) + 1;

    pub fn new(samp_rate: usize, center_freq: usize) -> AveragePsd<FFT_SIZE, AVG_WINDOW_SIZE> {
        let fft = FftPlanner::new().plan_fft_forward(FFT_SIZE);
        let psd_averages = vec!(NoSumSMA::<f32, f32, AVG_WINDOW_SIZE>::new();  FFT_SIZE);

        AveragePsd::<FFT_SIZE, AVG_WINDOW_SIZE>{samp_rate, center_freq, fft, psd_averages}
    }

    pub fn update(&mut self, samples: &mut [Complex<f32>; FFT_SIZE]) {
        let psd = self.calc_psd(samples).expect("unable to calculate PSD");
        // keep a running average of PSD values
        for (index, element) in psd.into_iter().enumerate() {
            self.psd_averages[index].add_sample(element);
        }
    }

    pub fn get_psd(&self) -> [f32; FFT_SIZE] {
        let mut avg_psd: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
        for (index, element) in self.psd_averages.iter().enumerate() {
            avg_psd[index] = element.get_average();
        }

        return avg_psd;
    }

    fn calc_psd(&self, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<[f32; FFT_SIZE]> {
        // https://pysdr.org/content/sampling.html#calculating-power-spectral-density

        // fft[0] -> DC (or samp_rate)
        // fft[len/2 + 1] -> nyquist f NOTE this assumes non-quadrature
        // fft[1] to fft[len/2] are positive frequencies, step is samp_rate/len(fft)
        // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/

        // both of these modify samples in-place
        self.hamming_window(samples).expect("failed to apply hamming window to samples");
        self.fft.process(samples);

        // "RustFFT does not normalize outputs. Callers must manually normalize
        // the results by scaling each element by 1/len().sqrt()"
        // https://docs.rs/rustfft/latest/rustfft/index.html#normalization
        let norm_factor: f32 = 1.0 / (FFT_SIZE as f32).sqrt();

        // TODO:reliability improve this
        assert_eq!(FFT_SIZE % 2, 0);

        // TODO:feature this should go from -samp_rate/2 to samp_rate/2 rather
        // than 0 - samp_rate, see https://pysdr.org/content/rtlsdr.html
        let psd: Vec<f32> = samples.into_iter().map(|fft_step| {
            (fft_step.norm() * norm_factor).log10() * 10.0 }).collect();

        let mut retval: [f32; FFT_SIZE] = psd.try_into().expect("failed convert vec->arr in calc_psd");
        Self::fftshift(&mut retval);

        Ok(retval)
    }

    fn hamming_window(&self, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<()> {
        // https://math.stackexchange.com/questions/248849/hamming-window-understanding-formula
        // TODO:optimization this really oughta be computed and cached inside
        // a struct for the fft - write the window as a const fn or something
        let mut window: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
        for n in 0..FFT_SIZE {
            let numerator = 2.0 * PI * (n as f32);
            window[n] = 0.53836 - 0.46164 * (numerator / (FFT_SIZE as f32 - 1.0)).cos();
        }

        // hamming formula is:
        // x = 0.53836 + 0.46164*cos( (2*pi*n) / (N - 1)) where N is fft_len

        for i in 0..FFT_SIZE {
            samples[i] = samples[i].scale(window[i]);
        }

        Ok(())
    }

    // https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html
    fn fftshift(fftbuf: &mut [f32; FFT_SIZE]) {

        // NOTE only works for even-length FFTs right now
        // TODO:reliability return better error, do compile time check, or
        // fix for odd len fftbuf
        assert_eq!(FFT_SIZE % 2, 0);

        // In Octave:
        // horzcat(f_ham(ceil((buflen-1)/2):-1:1), 0, f_ham(1:floor((buflen-1)/2)));
        let cutoff_point = ((FFT_SIZE-1) as f32)/2.0;
        fftbuf.rotate_left(cutoff_point.ceil() as usize);
    }

    pub fn get_freq_range(&self) -> [f32; FFT_SIZE] {
        // returns the frequency steps of the fftshift-ed FFT/PSD
        //let mut freq_range: [f32; (FFT_SIZE/2)+1] = [0.0; (FFT_SIZE/2)+1];
        let mut freq_range: [f32; FFT_SIZE] = [0.; FFT_SIZE];
        let mut i = 0;

        let delta_f: f32 = self.samp_rate as f32 / FFT_SIZE as f32;
        let start: f32 = self.center_freq as f32 - (self.samp_rate as f32 / 2.);
        while i < FFT_SIZE {
            freq_range[i] = start + (delta_f * i as f32);
            i += 1;
        }

        return freq_range;
    }
}


#[cfg(test)]
mod tests {
    use crate::AveragePsd;
    use num_complex::Complex;
    use wavegen;

    #[test]
    fn fftshift_even_len() {
        let mut startbuf: [f32; 8] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        // avg window of 1, FFT size of 8
        AveragePsd::<8, 1>::fftshift(&mut startbuf);
        let endbuf: [f32; 8] = [4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0, 3.0];
        assert_eq!(startbuf, endbuf);
    }

    #[test]
    #[should_panic(expected = "assertion `left == right` failed")]
    fn fftshift_odd_len() {
        let mut startbuf: [f32; 7] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // avg window of 1, FFT size of 8
        AveragePsd::<7,1>::fftshift(&mut startbuf);
        // should panic on the line above due to an assert() in fftshift
        let endbuf: [f32; 7] = [4.0, 5.0, 6.0, 0.0, 1.0, 2.0, 3.0];
        // this is not great but for now it'll work
        assert_eq!(startbuf, endbuf);
    }

    #[test]
    fn hamming_window_correct() {
        // avg window of 10, FFT size of 8
        const FFT_SIZE: usize = 8;
        let psd = AveragePsd::<FFT_SIZE, 10>::new(1, 1);

        let mut samples: [Complex<f32>; FFT_SIZE] = [Complex::new(1.0, 1.0); FFT_SIZE];
        psd.hamming_window(&mut samples).expect("failed to apply hamming window");

        let expected: [Complex<f32>; FFT_SIZE] = [Complex { re: 0.07672, im: 0.07672 }, Complex { re: 0.25053218, im: 0.25053218 }, Complex { re: 0.64108455, im: 0.64108455 }, Complex { re: 0.95428324, im: 0.95428324 }, Complex { re: 0.95428324, im: 0.95428324 }, Complex { re: 0.6410844, im: 0.6410844 }, Complex { re: 0.2505322, im: 0.2505322 }, Complex { re: 0.07672, im: 0.07672 }];

        assert_eq!(samples, expected);
    }

    #[test]
    fn all_zero_psd() {
        // pass all-zero samples into the code and get -inf decibels out
        // avg window of 10, FFT size of 1024
        const FFT_SIZE: usize = 1024;
        let psd = AveragePsd::<FFT_SIZE,10>::new(1, 1);

        let mut samples: [Complex<f32>; FFT_SIZE] = [Complex::new(0.0, 0.0); FFT_SIZE];

        let output = psd.calc_psd(&mut samples).expect("failed to calculate PSD");
        println!("{:?}", output);

        assert!(!output.contains(&f32::NAN));

        for elem in output.iter() {
            // log10(0) is -inf so -inf dB is correct, I think
            assert_eq!(*elem, std::f32::NEG_INFINITY);
        }
    }

    #[test]
    fn psd_of_sine_wave() {
        const FFT_SIZE: usize = 256;

        let samp_rate: f32 = 256.; // in Hz
        let center_freq: usize = 0;

        let mut psd = AveragePsd::<FFT_SIZE, 1>::new(samp_rate as usize, center_freq);

        let mut signal = wavegen::Waveform::<f32>::new(samp_rate);
        signal.add_component(wavegen::sine!(50.));

        let samples_re = signal.iter().take(FFT_SIZE).collect::<Vec<f32>>();
        let mut samples: [Complex<f32>; FFT_SIZE] = [Complex::new(0.0, 0.0); FFT_SIZE];

        for i in 0..samples_re.len() {
            samples[i] = Complex::new(samples_re[i], 0.0);
        }

        println!("Real samples:");
        println!("{:?}", samples_re);

        println!("Imaginary samples:");
        println!("{:?}", samples);

        psd.update(&mut samples);

        let elements = psd.get_psd();

        assert!(!elements.contains(&f32::NAN));

        println!("PSD:");
        println!("{:?}", elements);

        println!("PSD frequency steps:");
        println!("{:?}", psd.get_freq_range());

        assert_eq!(psd.get_freq_range().len(), FFT_SIZE);
        assert_eq!(psd.get_freq_range()[0], -128.);
        assert_eq!(psd.get_freq_range()[64], -64.);
        assert_eq!(psd.get_freq_range()[128], 0.);
        assert_eq!(psd.get_freq_range()[255], 127.);

        assert_eq!(elements.len(), FFT_SIZE);
        // make sure test frequency (50Hz) is present
        assert_eq!(psd.get_freq_range()[178], 50.);
        assert_eq!(elements[178], 1.6327056);

        // mirror freq at -50Hz should also be present
        assert_eq!(psd.get_freq_range()[78], -50.);
        assert_eq!(elements[78], 1.6327056);
    }

    // TODO:test test center_freq != 0
}
