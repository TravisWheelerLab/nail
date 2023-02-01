pub struct CloudSearchParams {
    pub target_start: usize,
    pub target_end: usize,
    pub profile_start: usize,
    pub profile_end: usize,
    pub gamma: usize,
    pub alpha: f32,
    pub beta: f32,
}

impl Default for CloudSearchParams {
    fn default() -> Self {
        CloudSearchParams {
            target_start: 0,
            target_end: 0,
            profile_start: 0,
            profile_end: 0,
            gamma: 5,
            alpha: 12.0,
            beta: 12.0,
        }
    }
}
