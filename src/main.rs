//#![feature(generic_const_exprs)]

mod averagepsd;
mod monitor;
mod client;

use std::thread;
use std::net::TcpListener;
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

fn main() {
    // TODO take in config
    // probably YAML or something defining the radios + bands
    let dev_filter: String = "driver=rtlsdr".to_string();
    let center_freq: usize = 86_200_000; // where to tune the SDR

    // enumerate the available SDR devices
    let devs = soapysdr::enumerate(&dev_filter[..]).expect("Error listing devices");
    let dev_args = match devs.len() {
        0 => {
            eprintln!("no matching SDR devices found");
            std::process::exit(1);
        }
        1 => devs.into_iter().next().unwrap(),
        n => {
            eprintln!("{} devices found. Choose from:",n);
            for dev in devs {
                eprintln!("\t'{}'", dev);
            }
            std::process::exit(1);
        }
    };

    // set up channel - cloned for each monitor
    let (mon_send, mon_recv) = crossbeam_channel::unbounded();

    // create monitor(s)
    // TODO only works for one monitor right now
    let dev = soapysdr::Device::new(dev_args).expect("Error opening device");
    let mut mon = Monitor::<FFT_SIZE, AVG_WINDOW_SIZE>::new(dev, mon_send.clone(), mon_recv.clone(), SAMP_RATE, center_freq);
    // start monitor(s)
    // TODO keep the handles, join() on shutdown
    let _handle = thread::spawn(move || mon.start());

    // start TCP server

    let listener = TcpListener::bind("127.0.0.1:8080")
        .expect("should be able to bind to a local port");

    println!("TCP listener started on 127.0.0.1:8080"); // TODO

    let (tcp_send, tcp_rcv) = crossbeam_channel::unbounded();
    let _handle = thread::spawn(move || serve_tcp(listener, tcp_send));

    let mut clients: Vec<crossbeam_channel::Sender<Message>> = Vec::new();
    loop {
        // get any new clients
        match tcp_rcv.try_recv() {
            Ok(mut new_client) => {
                clients.push(new_client.get_send_channel());
                let _handle = thread::spawn(move || new_client.start());
                println!("Client successfully subscribed.");
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
                //println!("{:?}", serde_json::to_value(&data).unwrap());
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
    println!("TCP up and running...");
    for stream in listener.incoming() {
        match stream {
            Ok(stream) => {
                let addr = stream.peer_addr().expect("should have peer address in TCP stream");
                println!("new client connected: {addr:?}");

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
