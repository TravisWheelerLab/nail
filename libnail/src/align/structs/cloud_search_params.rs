#[derive(Clone)]
pub struct CloudSearchParams {
    pub gamma: usize,
    pub alpha: f32,
    pub beta: f32,
}

impl Default for CloudSearchParams {
    fn default() -> Self {
        CloudSearchParams {
            gamma: 5,
            alpha: 12.0,
            beta: 20.0,
        }
    }
}
