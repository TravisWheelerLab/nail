use anyhow::Result;
use serde::Serialize;
use std::io::Write;

#[derive(Serialize)]
pub struct SodaAnnotation {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub row: usize,
}

pub trait SodaJson {
    fn soda_json(&self, out: &mut impl Write) -> Result<()>;
}
