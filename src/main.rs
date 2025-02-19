//#![feature(generic_const_exprs)]

mod averagepsd;
mod monitor;

use nbfm_monitor_ui::messages::Message;
use crate::monitor::Monitor;

const FFT_SIZE: usize = 512;
// TODO figure out how to vary this (or hardcode per-SDR-model? idk)
const SAMP_RATE: usize = 2_560_000; // 2.56 MHz (for rtlsdr)
// window of the running average
const AVG_WINDOW_SIZE: usize = (SAMP_RATE / FFT_SIZE) * 2; // 2 seconds of averaging

fn main() {
    // TODO take in config
    let dev_filter: String = "driver=rtlsdr".to_string();
    let center_freq: usize = 86_200_000; // where to tune the SDR

    // enumerate the available SDR devices
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

    // set up channel - cloned for each monitor
    let (send, recv) = crossbeam_channel::unbounded();

    // create monitor(s)
    let dev = soapysdr::Device::new(dev_args).expect("Error opening device");
    let (send2, recv2) = (send.clone(), recv.clone());

    let mut mon = Monitor::<FFT_SIZE, AVG_WINDOW_SIZE>::new(dev, send2, recv2, SAMP_RATE, center_freq);
    // start TCP server
    // start monitor(s)
    // TODO in threads
    mon.start();
    // route data to TCP clients

    // eventually, stop
    send.send(Message::Stop{}).expect("Should be able to send a message to stop threads");
}


#[cfg(test)]
mod tests {
    #[test]
    fn runs_at_all() {
        assert_eq!(1,1);
    }

}
