//#![feature(generic_const_exprs)]

mod averagepsd;
mod monitor;
mod client;

use config::Config;
use chrono::Local;
use colog;
use colog::format::CologStyle;
use log::{Level, LevelFilter, error, warn, info, debug, trace};
use colored::Colorize;
use std::thread;
use std::net::TcpListener;
use std::collections::HashMap;
use crossbeam_channel::{TryRecvError, TrySendError};

use nbfm_monitor_ui::messages::Message;
use crate::monitor::Monitor;
use crate::client::Client;

const FFT_SIZE: usize = 512;
// TODO figure out how to vary this (or hardcode per-SDR-model? idk)
// could use the bandwidth param from soapysdr
const SAMP_RATE: usize = 2_560_000; // 2.56 MHz (for rtlsdr)
// window of the running average
const AVG_WINDOW_SIZE: usize = (SAMP_RATE / FFT_SIZE) * 2; // 2 seconds of averaging

// configure logging output style
pub struct CustomPrefixToken;
impl CologStyle for CustomPrefixToken {
    fn prefix_token(&self, level: &Level) -> String {
        let now = Local::now().format("%H:%M:%S");
        format!(
            "{}{} {}{}",
            "[".blue().bold(),
            self.level_color(level, self.level_token(level)),
            now,
            "]".blue().bold(),
        )
    }

    fn level_token(&self, level: &Level) -> &str {
        match *level {
            Level::Error => "ERR",
            Level::Warn => "WRN",
            Level::Info => "INF",
            Level::Debug => "DBG",
            Level::Trace => "TRC",
        }
    }
}

fn main() {
    // set up logging
    let mut builder = colog::default_builder();
    builder.format(colog::formatter(CustomPrefixToken));
    builder.filter(None, LevelFilter::Trace);
    builder.init();

    // load config, set defaults
    // probably YAML or something defining the radios + bands
    let settings = Config::builder()
        .set_default("dev_filter", "driver=rtlsdr").unwrap()
        .set_default("center_freq", 86_200_000).unwrap()
        .add_source(config::File::with_name("settings"))
        .build()
        .unwrap();

    info!("Started with settings:\n\t{:?}",
        settings.clone().try_deserialize::<HashMap<String, String>>()
        .unwrap()
    );

    // enumerate the available SDR devices
    let dev_filter: String = settings.get("dev_filter").unwrap();
    let devs = soapysdr::enumerate(&dev_filter[..]).expect("Error listing devices");
    let dev_args = match devs.len() {
        0 => {
            error!("no matching SDR devices found");
            std::process::exit(1);
        }
        1 => devs.into_iter().next().unwrap(),
        n => {
            let mut errlog = String::new().to_owned();
            for dev in devs {
                errlog = format!("{errlog}\t'{dev}'");
            }
            error!("{} devices found. Choose from: {}", n, errlog);
            std::process::exit(1);
        }
    };

    // set up channel - cloned for each monitor
    let (mon_send, mon_recv) = crossbeam_channel::unbounded();

    // create monitor(s)
    // TODO only works for one monitor right now
    let dev = soapysdr::Device::new(dev_args).expect("Error opening device");
    let center_freq: usize = settings.get("center_freq").unwrap();
    let mut mon = Monitor::<FFT_SIZE, AVG_WINDOW_SIZE>::new(
        dev,
        mon_send.clone(),
        mon_recv.clone(),
        SAMP_RATE,
        center_freq,
    );

    // start monitor(s)
    // TODO keep the handles, join() on shutdown
    let _handle = thread::spawn(move || mon.start());

    // start TCP server

    let listener = TcpListener::bind("127.0.0.1:8080")
        .expect("should be able to bind to a local port");

    info!("TCP listener started on 127.0.0.1:8080"); // TODO

    let (tcp_send, tcp_rcv) = crossbeam_channel::unbounded();
    let _handle = thread::spawn(move || serve_tcp(listener, tcp_send));

    let mut clients: Vec<crossbeam_channel::Sender<Message>> = Vec::new();
    loop {
        // get any new clients
        match tcp_rcv.try_recv() {
            Ok(mut new_client) => {
                clients.push(new_client.get_send_channel());
                let _handle = thread::spawn(move || new_client.start());
                info!("Client successfully subscribed.");
            }
            Err(TryRecvError::Empty) => {},
            Err(e) => panic!("{e}"),
        }

        // get new frames from the monitors
        let msg = mon_recv.recv().expect("receive message from monitors");
        match msg {
            Message::Frame{ .. } => {
                // we got a frame, forward it to the clients
                clients.retain(|client| {
                    // return true to keep element
                    // TODO get rid of clone here? (or at least change to copy)
                    return match client.try_send(msg.clone()) {
                        Ok(_) => true,
                        Err(TrySendError::Disconnected(_)) => false,
                        Err(e) => panic!("message send panic: {e:?}"),
                    }
                });
                //debug!("{:?}", serde_json::to_value(&data).unwrap());
            }
            Message::Connect{ .. } => {}
            Message::ConnectResult{ .. } => {}
            Message::Error{ .. } => {}
            Message::Stop {} => {}
        }
    }

    // eventually, stop
    // TODO catch ctrl-c and:
    //   - send STOP to each monitor
    //   - send STOP to TCP listener thread
    //   - wait for each thread to join()
    //
    // send.send(Message::Stop{}).expect("Should be able to send a message to stop threads");
    // for handle in handles: handle.join().unwrap();
}

fn serve_tcp(listener: TcpListener, sender: crossbeam_channel::Sender<Client>) {
    info!("TCP up and running...");
    for stream in listener.incoming() {
        match stream {
            Ok(stream) => {
                let addr = stream.peer_addr().expect("should have peer address in TCP stream");
                info!("new client connected: {addr:?}");

                let client = Client::new(stream);
                sender.send(client).expect("send new client to main thread");
            }
            Err(e) => match e.kind() {
                _ => panic!("IO error while TCP connecting: {e}"),
            }
        }
    }
    // cleanup here
}

#[cfg(test)]
mod tests {
    #[test]
    fn runs_at_all() {
        assert_eq!(1,1);
    }
}
