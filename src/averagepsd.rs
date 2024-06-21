use std::f32::consts::PI;
use rustfft::{FftPlanner, Fft};
use simple_moving_average::{SMA, NoSumSMA};
use std::sync::Arc;
use num_complex::Complex;
use std::io;

// keeps a running average of the power spectral density of a frequency range
pub struct AveragePsd<const AVG_WINDOW_SIZE: usize, const FFT_SIZE: usize> {
    samp_rate: usize,
    center_freq: usize,
    fft: Arc<dyn Fft<f32>>,
    psd_averages: Vec<NoSumSMA::<f32, f32, AVG_WINDOW_SIZE>>,
}

impl<const AVG_WINDOW_SIZE: usize, const FFT_SIZE: usize> AveragePsd<AVG_WINDOW_SIZE, FFT_SIZE> {
    pub fn new(samp_rate: usize, center_freq: usize) -> AveragePsd<AVG_WINDOW_SIZE, FFT_SIZE> {
        let fft = FftPlanner::new().plan_fft_forward(FFT_SIZE);
        let psd_averages = vec!(NoSumSMA::<f32, f32, AVG_WINDOW_SIZE>::new(); FFT_SIZE);

        AveragePsd{samp_rate, center_freq, fft, psd_averages}
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

    // TODO don't use io::Error here, I guess, probably
    fn calc_psd(&self, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<[f32; FFT_SIZE]> {
        // https://pysdr.org/content/sampling.html#calculating-power-spectral-density

        // fft[0] -> DC
        // fft[len/2 + 1] -> nyquist f
        // fft[1] to fft[len/2] are positive frequencies, step is samp_rate/len(fft)
        // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/

        // both of these modify samples in-place
        self.hamming_window(samples).expect("failed to apply hamming window to samples");

        self.fft.process(samples);

        // output will go here (f32 instead of complex like the input)
        let mut psd: [f32; FFT_SIZE] = [0.0; FFT_SIZE];

        // "RustFFT does not normalize outputs. Callers must manually normalize
        // the results by scaling each element by 1/len().sqrt()"
        // https://docs.rs/rustfft/latest/rustfft/index.html#normalization
        let norm_factor:f32 = 1.0 / (FFT_SIZE as f32).sqrt();
        for i in 0..FFT_SIZE {
            // normalize FFT output
            psd[i] = samples[i].norm() * norm_factor;
            psd[i] = 10.0 * psd[i].log10(); // convert to dB
        }

        Self::fftshift(&mut psd);
        // psd now contains our completed PSD calculation

        Ok(psd)
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
        // returns the frequency steps of the FFT/PSD
        let mut freq_range: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
        let mut i = 0;

        let start_f: f32 = (self.samp_rate as f32 / -2.0) + self.center_freq as f32;
        let delta_f: f32 = self.samp_rate as f32 / FFT_SIZE as f32;
        while i < FFT_SIZE {
            freq_range[i] = start_f + (delta_f * i as f32);
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
        AveragePsd::<1, 8>::fftshift(&mut startbuf);
        let endbuf: [f32; 8] = [4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0, 3.0];
        assert_eq!(startbuf, endbuf);
    }

    #[test]
    #[should_panic(expected = "assertion `left == right` failed")]
    fn fftshift_odd_len() {
        let mut startbuf: [f32; 7] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // avg window of 1, FFT size of 8
        AveragePsd::<1, 7>::fftshift(&mut startbuf);
        // should panic on the line above due to an assert() in fftshift
        let endbuf: [f32; 7] = [4.0, 5.0, 6.0, 0.0, 1.0, 2.0, 3.0];
        // this is not great but for now it'll work
        assert_eq!(startbuf, endbuf);
    }

    #[test]
    fn hamming_window_correct() {
        // avg window of 10, FFT size of 8
        const FFT_SIZE: usize = 8;
        let psd = AveragePsd::<10, FFT_SIZE>::new(1, 1);

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
        let psd = AveragePsd::<10, FFT_SIZE>::new(1, 1);

        let mut samples: [Complex<f32>; FFT_SIZE] = [Complex::new(0.0, 0.0); FFT_SIZE];

        let output = psd.calc_psd(&mut samples).expect("failed to calculate PSD");
        println!("{:?}", output);

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

        let mut psd = AveragePsd::<1, FFT_SIZE>::new(samp_rate as usize, center_freq);


        let mut signal = wavegen::Waveform::<f32>::new(samp_rate);
        signal.add_component(wavegen::sine!(50.));

        let samples_re = signal.iter().take(FFT_SIZE).collect::<Vec<f32>>();
        let mut samples: [Complex<f32>; FFT_SIZE] = [Complex::new(0.0, 0.0); FFT_SIZE];

        for i in 0..samples_re.len() {
            samples[i] = Complex::new(samples_re[i], 0.0);
        }

        println!("{:?}", samples_re);
        println!("{:?}", samples);

        psd.update(&mut samples);

        let elements = psd.get_psd();
        println!("{:?}", elements);
        assert_eq!(elements.len(), FFT_SIZE);
        // symmetric peaks
        assert!((elements[79] - elements[179]).abs() < 0.0001);

        // peaks should be ~= -2.67 dB
        assert!(elements[79] > 2.6);
        assert!(elements[179] < 2.7);

        assert!(elements[79] > 2.6);
        assert!(elements[179] < 2.7);

        println!("{:?}", psd.get_freq_range());
        assert_eq!(psd.get_freq_range().len(), FFT_SIZE);
        // TODO:bugfix I think this is actually wrong, should start at -127 I
        // think? (missing one of the two removed values of the FFT maybe?)
        assert_eq!(psd.get_freq_range()[0], -128.0);

        // bins 128 - 255 are 0-127hz
        assert_eq!(psd.get_freq_range()[128], center_freq as f32);
        assert_eq!(psd.get_freq_range()[255], 127.0);

    }
}
