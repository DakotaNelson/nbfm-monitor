use serde::Deserialize;

// expected structure of the config file
#[derive(Debug, Deserialize)]
#[allow(unused)]
pub struct AppConfig {
    pub radios: Vec<RadioConfig>,
    pub tcp: TcpConfig,
}


#[derive(Debug, Deserialize)]
#[allow(unused)]
pub struct RadioConfig {
    pub filter: String,
    pub freq: f64, // in MHz - DO NOT forget to convert!
}

#[derive(Debug, Deserialize)]
#[allow(unused)]
pub struct TcpConfig {
    pub ip: String,
    pub port: u16,
}
