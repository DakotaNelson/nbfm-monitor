use std::f32::consts::PI;
use rustfft::{FftPlanner, Fft};
use rustfft::num_traits::Pow;
use simple_moving_average::{SMA, NoSumSMA};
use std::sync::Arc;
use num_complex::Complex;
use std::io;

// keeps a running average of the power spectral density of a frequency range
pub struct AveragePsd<const AVG_WINDOW_SIZE: usize, const FFT_SIZE: usize> {
    samp_rate: usize,
    CENTER_FREQ: usize,
    //FFT_SIZE: usize,
    //AVG_WINDOW_SIZE: usize,
    fft: Arc<dyn Fft<f32>>,
    pub fft_averages: Vec<NoSumSMA::<f32, f32, {AVG_WINDOW_SIZE}>>,
    // TODO write a getter for this
}

impl<const AVG_WINDOW_SIZE: usize, const FFT_SIZE: usize> AveragePsd<AVG_WINDOW_SIZE, FFT_SIZE> {
    pub fn new(samp_rate: usize, CENTER_FREQ: usize) -> AveragePsd<AVG_WINDOW_SIZE, FFT_SIZE> {
        let fft = FftPlanner::new().plan_fft_forward(FFT_SIZE);
        let fft_averages = vec!(NoSumSMA::<f32, f32, { AVG_WINDOW_SIZE }>::new(); FFT_SIZE);

        AveragePsd{samp_rate, CENTER_FREQ,
                   fft, fft_averages}
    }

    pub fn update(&mut self, samples: &mut [Complex<f32>; FFT_SIZE]) {
        let psd = self.calc_psd(samples).expect("unable to calculate PSD");
        // keep a running average of PSD values
        for (index, element) in psd.into_iter().enumerate() {
            self.fft_averages[index].add_sample(element);
        }
    }

    // TODO don't use io::Error here, I guess, probably
    fn calc_psd(&self, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<[f32; FFT_SIZE]> {
        // https://pysdr.org/content/sampling.html#calculating-power-spectral-density

        // fft[0] -> DC
        // fft[len/2 + 1] -> nyquist f
        // 2 - len/2 are positive frequencies, each step is samp_rate/len(fft)
        // https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/

        // both of these modify samples in-place
        self.hamming_window(samples).expect("failed to apply hamming window to samples");
        self.fft.process(samples);

        let mut psd: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
        for i in 0..FFT_SIZE {
            let sample: f32 = samples[i].norm().pow(2) / (FFT_SIZE * self.samp_rate) as f32;
            psd[i] = 10.0 * sample.log10(); // convert to dB
        }

        Self::fftshift(&mut psd.to_vec());
        // psd now contains our completed PSD calculation

        Ok(psd)
    }

    fn hamming_window(&self, samples: &mut [Complex<f32>; FFT_SIZE]) -> io::Result<()> {
        // https://math.stackexchange.com/questions/248849/hamming-window-understanding-formula
        // TODO this really oughta be computed and cached inside a struct for the fft
        // TODO - write the window as a const fn or something
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

    pub fn get_freq_range(&self) -> [f32; FFT_SIZE] {
        // returns the frequency steps of the FFT/PSD
        let mut freq_range: [f32; FFT_SIZE] = [0.0; FFT_SIZE];
        let mut i = 0;

        while i < FFT_SIZE {
            freq_range[i] = self.CENTER_FREQ as f32 + (self.samp_rate as f32/-2.0) + ((self.samp_rate/FFT_SIZE) * i) as f32;
            i += 1;
        }

        return freq_range;
    }
}


#[cfg(test)]
mod tests {
    use crate::AveragePsd;

    #[test]
    fn fftshift_works() {
        let mut startbuf: Vec<f32> = vec![0.0, 1.0, 2.0, 3.0, -4.0, -3.0, -2.0, -1.0];
        // avg window of 10, FFT size of 8
        AveragePsd::<10, 8>::fftshift(&mut startbuf);
        let endbuf: Vec<f32> = vec![-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        assert_eq!(startbuf, endbuf);
    }
}
