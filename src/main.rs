//#![feature(generic_const_exprs)]

mod averagepsd;
mod monitor;

//use nbfm_monitor_ui::messages::Message;
use crate::monitor::Monitor;

const FFT_SIZE: usize = 512;
// TODO figure out how to vary this (or hardcode per-SDR-model? idk)
const SAMP_RATE: usize = 2_560_000; // 2.56 MHz (for rtlsdr)
// window of the running average
const AVG_WINDOW_SIZE: usize = (SAMP_RATE / FFT_SIZE) * 1; // 2 seconds of averaging
                                                       // (only compute every *other*
                                                       // FFT right now)

fn main() {
    // take in config
    let dev_filter: String = "driver=rtlsdr".to_string();
    let center_freq: usize = 86_200_000; // where to tune the SDR
    // create monitor(s)
    let mut mon = Monitor::<FFT_SIZE, AVG_WINDOW_SIZE>::new(SAMP_RATE, center_freq, dev_filter);
    // start TCP server
    // start monitor(s)
    mon.start();
    // route data to clients
}


#[cfg(test)]
mod tests {
    #[test]
    fn runs_at_all() {
        assert_eq!(1,1);
    }

}
