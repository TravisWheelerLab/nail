use crate::util::PrintMe;

pub struct Bound {
    pub target_idx: usize,
    pub profile_idx: usize,
}

impl PrintMe for Bound {
    fn print(&self) {
        println!("{}-{}", self.target_idx, self.profile_idx);
    }
}
