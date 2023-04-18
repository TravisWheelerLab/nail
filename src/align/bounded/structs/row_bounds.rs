use crate::align::bounded::structs::bound::CloudBoundGroup;
use anyhow::Result;
use std::fmt::{Debug, Formatter};
use std::io::Write;

pub struct RowBounds {
    pub target_start: usize,
    pub target_end: usize,
    pub row_capacity: usize,
    pub left_row_bounds: Vec<usize>,
    pub right_row_bounds: Vec<usize>,
}

impl Default for RowBounds {
    fn default() -> Self {
        Self {
            target_start: 1,
            target_end: 1,
            row_capacity: 1,
            left_row_bounds: vec![1, 1],
            right_row_bounds: vec![1, 1],
        }
    }
}

impl RowBounds {
    pub fn new(cloud_bounds: &CloudBoundGroup) -> Self {
        let mut params = Self::default();
        params.reuse(cloud_bounds);
        params
    }

    pub fn reuse(&mut self, cloud_bounds: &CloudBoundGroup) {
        let first_bound = cloud_bounds.get_first();
        let last_bound = cloud_bounds.get_last();

        self.target_start = first_bound.right_target_idx;
        self.target_end = last_bound.left_target_idx;

        let num_rows = cloud_bounds.target_length + 1;

        if num_rows > self.row_capacity {
            self.left_row_bounds.resize(num_rows, usize::MAX);
            self.right_row_bounds.resize(num_rows, usize::MIN);
            self.row_capacity = num_rows;
        }

        for bound in &cloud_bounds.bounds {
            // TODO: this is a band-aid fix!
            if bound.was_pruned() {
                continue;
            }

            self.left_row_bounds[bound.left_target_idx] =
                self.left_row_bounds[bound.left_target_idx].min(bound.left_profile_idx);

            self.left_row_bounds[bound.right_target_idx] =
                self.left_row_bounds[bound.right_target_idx].min(bound.right_profile_idx);

            self.right_row_bounds[bound.left_target_idx] =
                self.right_row_bounds[bound.left_target_idx].max(bound.left_profile_idx);

            self.right_row_bounds[bound.right_target_idx] =
                self.right_row_bounds[bound.right_target_idx].max(bound.right_profile_idx);
        }
    }

    pub fn valid(&self) -> bool {
        for row_idx in self.target_start..=self.target_end {
            if self.left_row_bounds[row_idx] > self.right_row_bounds[row_idx] {
                return false;
            }
        }
        true
    }

    pub fn dump(&self, out: &mut impl Write) -> Result<()> {
        for row_idx in self.target_start..=self.target_end {
            writeln!(
                out,
                "{}: {}-{}",
                row_idx, self.left_row_bounds[row_idx], self.right_row_bounds[row_idx]
            )?;
        }
        Ok(())
    }
}

impl Debug for RowBounds {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "rows: {}-{}", self.target_start, self.target_end)?;
        for row_idx in self.target_start..=self.target_end {
            writeln!(
                f,
                "{}: {}-{}",
                row_idx, self.left_row_bounds[row_idx], self.right_row_bounds[row_idx]
            )?;
        }
        Ok(())
    }
}
