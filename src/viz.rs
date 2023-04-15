use serde::Serialize;

#[derive(Serialize)]
pub struct SodaAnnotation {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub row: usize,
}
