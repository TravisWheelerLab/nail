use crate::util::PrintMe;
use crate::viz::{SodaAnnotation, SodaJson};
use anyhow::Result;
use serde::Serialize;
use std::io::Write;
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

#[derive(Clone)]
pub struct CloudBound {
    pub left_target_idx: usize,
    pub left_profile_idx: usize,
    pub right_target_idx: usize,
    pub right_profile_idx: usize,
}

impl Default for CloudBound {
    fn default() -> Self {
        // set the default to a simple invalid bound
        // (the left bound is on the right, and vice versa)
        CloudBound {
            left_target_idx: 0,
            left_profile_idx: 1,
            right_target_idx: 1,
            right_profile_idx: 0,
        }
    }
}

impl CloudBound {
    pub fn was_pruned(&self) -> bool {
        self.left_target_idx < self.right_target_idx
    }

    pub fn anti_diagonal_idx(&self) -> usize {
        self.left_target_idx + self.left_profile_idx
    }

    pub fn anti_diagonal_cell_zip(&self) -> Zip<Rev<RangeInclusive<usize>>, RangeInclusive<usize>> {
        let target_range = (self.right_target_idx..=self.left_target_idx).rev();
        let profile_range = self.left_profile_idx..=self.right_profile_idx;
        target_range.zip(profile_range)
    }
}

pub struct CloudBoundGroup {
    pub bounds: Vec<CloudBound>,
    pub size: usize,
    pub min_anti_diagonal_idx: usize,
    pub max_anti_diagonal_idx: usize,
}

impl CloudBoundGroup {
    pub fn new(size: usize) -> Self {
        Self {
            bounds: vec![CloudBound::default(); size],
            size,
            min_anti_diagonal_idx: size,
            max_anti_diagonal_idx: 0,
        }
    }

    pub fn resize(&mut self, new_size: usize) {
        self.bounds.resize(new_size, CloudBound::default());
        self.size = new_size;
    }

    pub fn reuse(&mut self, size: usize) {
        if size > self.size {
            self.resize(size);
        }

        self.min_anti_diagonal_idx = size;
        self.max_anti_diagonal_idx = 0;

        for bound in self.bounds.iter_mut() {
            bound.left_target_idx = 0;
            bound.left_profile_idx = 1;
            bound.right_target_idx = 1;
            bound.right_profile_idx = 0;
        }
    }

    pub fn set(
        &mut self,
        anti_diagonal_idx: usize,
        left_target_idx: usize,
        left_profile_idx: usize,
        right_target_idx: usize,
        right_profile_idx: usize,
    ) {
        self.min_anti_diagonal_idx = self.min_anti_diagonal_idx.min(anti_diagonal_idx);
        self.max_anti_diagonal_idx = self.max_anti_diagonal_idx.max(anti_diagonal_idx);
        let bound = &mut self.bounds[anti_diagonal_idx];
        bound.left_target_idx = left_target_idx;
        bound.left_profile_idx = left_profile_idx;
        bound.right_target_idx = right_target_idx;
        bound.right_profile_idx = right_profile_idx;
    }

    pub fn get(&self, idx: usize) -> &CloudBound {
        &self.bounds[idx]
    }

    pub fn get_mut(&mut self, idx: usize) -> &mut CloudBound {
        &mut self.bounds[idx]
    }
}

pub fn join_bounds(
    forward_bounds: &mut CloudBoundGroup,
    backward_bounds: &CloudBoundGroup,
) -> Result<()> {
    let start_idx = forward_bounds
        .min_anti_diagonal_idx
        .min(backward_bounds.min_anti_diagonal_idx);

    let end_idx = forward_bounds
        .max_anti_diagonal_idx
        .max(backward_bounds.max_anti_diagonal_idx);

    forward_bounds.min_anti_diagonal_idx = start_idx;
    forward_bounds.max_anti_diagonal_idx = end_idx;
    let forward_slice = &mut forward_bounds.bounds[start_idx..=end_idx];
    let backward_slice = &backward_bounds.bounds[start_idx..=end_idx];

    for (forward_bound, backward_bound) in forward_slice.iter_mut().zip(backward_slice) {
        if forward_bound.was_pruned() {
            // if there's no valid forward bound, just take the backward bound
            forward_bound.left_target_idx = backward_bound.left_target_idx;
            forward_bound.left_profile_idx = backward_bound.left_profile_idx;
            forward_bound.right_target_idx = backward_bound.right_target_idx;
            forward_bound.right_profile_idx = backward_bound.right_profile_idx;
        } else if backward_bound.was_pruned() {
            // if there's no valid backward bound, we can do nothing since we
            // are consuming the forward bounds
            continue;
        }
        // otherwise we have two valid bounds and we can compare them
        forward_bound.left_target_idx = forward_bound
            .left_target_idx
            .max(backward_bound.left_target_idx);

        forward_bound.left_profile_idx = forward_bound
            .left_profile_idx
            .min(backward_bound.left_profile_idx);

        forward_bound.right_target_idx = forward_bound
            .right_target_idx
            .min(backward_bound.right_target_idx);

        forward_bound.right_profile_idx = forward_bound
            .right_profile_idx
            .max(backward_bound.right_profile_idx);
    }

    Ok(())
}

impl PrintMe for CloudBound {
    fn print(&self) {
        println!(
            "{},{} : {},{}",
            self.left_profile_idx,
            self.left_target_idx,
            self.right_profile_idx,
            self.right_target_idx
        );
    }
}

impl PrintMe for Vec<CloudBound> {
    fn print(&self) {
        for b in self {
            b.print();
        }
    }
}

#[derive(Serialize)]
pub struct CloudBoundAnnotations {
    left: Vec<SodaAnnotation>,
    right: Vec<SodaAnnotation>,
}

impl SodaJson for CloudBoundGroup {
    fn soda_json(&self, out: &mut impl Write) -> Result<()> {
        let left: Vec<SodaAnnotation> = self.bounds
            [self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx]
            .iter()
            .map(|b| SodaAnnotation {
                id: format!("left-{}", b.left_target_idx + b.left_profile_idx),
                start: b.left_profile_idx,
                end: b.left_profile_idx + 1,
                row: b.left_target_idx,
            })
            .collect();

        let right: Vec<SodaAnnotation> = self.bounds
            [self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx]
            .iter()
            .map(|b| SodaAnnotation {
                id: format!("right-{}", b.right_target_idx + b.right_profile_idx),
                start: b.right_profile_idx,
                end: b.right_profile_idx + 1,
                row: b.right_target_idx,
            })
            .collect();

        let cloud_ann = CloudBoundAnnotations { left, right };
        let serialized = serde_json::to_string(&cloud_ann).unwrap();
        write!(out, "{serialized}")?;
        Ok(())
    }
}
