use serde::Deserialize;

// expected structure of the config file
#[derive(Debug, Deserialize)]
#[allow(unused)]
pub struct AppConfig {
    pub radios: Vec<RadioConfig>,
}


#[derive(Debug, Deserialize)]
#[allow(unused)]
pub struct RadioConfig {
    pub filter: String,
    pub freq: f64, // in MHz - DO NOT forget to convert!
}
