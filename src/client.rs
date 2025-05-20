// a connected TCP client that receives data from the monitoring system

use std::net::{TcpStream, SocketAddr};
use std::io::Write;
use log::{info};

use nbfm_monitor_ui::messages::Message;

pub struct Client {
    sender: crossbeam_channel::Sender<Message>,
    receiver: crossbeam_channel::Receiver<Message>,
    addr: SocketAddr,
    stream: TcpStream,
}

impl Client {
    pub fn new(stream: TcpStream) -> Client {
        let (send, rcv) = crossbeam_channel::unbounded();
        let peer_address = stream.peer_addr().expect("should be able to get peer addr from TCP stream");
        return Client {
            sender: send,
            receiver: rcv,
            addr: peer_address,
            stream: stream,
        }
    }

    pub fn get_send_channel(&self) -> crossbeam_channel::Sender<Message> {
        return self.sender.clone();
    }

    pub fn start(&mut self) {
        info!("Client connected");
        loop {
            let msg = self.receiver.recv().expect("should get a message");

            let mut writer = std::io::BufWriter::new(&mut self.stream);
            serde_json::to_writer(&mut writer, &msg)
                .expect("should be able to serialize frame message");

            // without flush, ~5 frames are batched to send, causing lag
            // maybe set_nodelay instead?
            // stream.set_nodelay(true).expect("set_nodelay call failed");
            match writer.flush() {
                Ok(_) => {},
                Err(e) => match e.kind() {
                    std::io::ErrorKind::BrokenPipe => {
                        info!("Client {0} disconnected (broken pipe)", self.addr);
                        break;
                    }
                    _ => panic!("failed to flush TCP socket: {e}"),
                }
            };
        }
    }
}
