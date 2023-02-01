use crate::util::PrintMe;
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

#[derive(Clone)]
pub struct CloudBound {
    pub left_target_idx: usize,
    pub left_profile_idx: usize,
    pub right_target_idx: usize,
    pub right_profile_idx: usize,
}

impl CloudBound {
    pub fn was_pruned(&self) -> bool {
        self.right_profile_idx < self.left_profile_idx
    }

    pub fn anti_diagonal(&self) -> Zip<Rev<RangeInclusive<usize>>, RangeInclusive<usize>> {
        let target_range = (self.right_target_idx..=self.left_target_idx).rev();
        let profile_range = self.left_profile_idx..=self.right_profile_idx;
        target_range.zip(profile_range)
    }
}

impl PrintMe for CloudBound {
    fn print(&self) {
        println!(
            "{},{} : {},{}",
            self.left_target_idx,
            self.left_profile_idx,
            self.right_target_idx,
            self.right_profile_idx
        );
    }
}
