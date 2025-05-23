use crate::align::structs::anti_diagonal_bounds::Cloud;
use crate::util::{MaxAssign, MinAssign};
use anyhow::Result;
use std::fmt::{Debug, Formatter};
use std::io::Write;

#[derive(Clone)]
pub struct RowBounds {
    pub seq_len: usize,
    pub prf_len: usize,
    pub seq_start: usize,
    pub seq_end: usize,
    pub row_capacity: usize,
    pub left_row_bounds: Vec<usize>,
    pub right_row_bounds: Vec<usize>,
    pub num_cells: usize,
}

impl Default for RowBounds {
    fn default() -> Self {
        Self {
            seq_len: 0,
            prf_len: 0,
            seq_start: 0,
            seq_end: 0,
            row_capacity: 1,
            left_row_bounds: vec![usize::MAX],
            right_row_bounds: vec![0],
            num_cells: 0,
        }
    }
}

impl RowBounds {
    pub fn new(target_end: usize) -> Self {
        let mut bounds = Self::default();
        bounds.reuse(target_end);
        bounds
    }

    pub fn reuse(&mut self, target_end: usize) {
        for row_idx in self.seq_start..=self.seq_end {
            self.left_row_bounds[row_idx] = usize::MAX;
            self.right_row_bounds[row_idx] = 0;
        }

        let num_rows = target_end + 1;

        self.left_row_bounds.resize(num_rows, usize::MAX);
        self.right_row_bounds.resize(num_rows, 0);
        self.row_capacity = num_rows;
        self.num_cells = 0;
    }

    pub fn fill_from_cloud(&mut self, cloud: &Cloud) {
        self.seq_len = cloud.seq_len;
        self.prf_len = cloud.prf_len;

        self.seq_start = cloud.first().0.seq_idx;
        self.seq_end = cloud.last().1.seq_idx;

        for bound in cloud.iter() {
            self.left_row_bounds[bound.1.seq_idx].min_assign(bound.1.prf_idx);
            self.left_row_bounds[bound.0.seq_idx].min_assign(bound.0.prf_idx);
            self.right_row_bounds[bound.1.seq_idx].max_assign(bound.1.prf_idx);
            self.right_row_bounds[bound.0.seq_idx].max_assign(bound.0.prf_idx);
            self.num_cells += bound.len()
        }

        // if the cloud bounds were set up properly, we should have no 0's
        (self.seq_start..=self.seq_end).for_each(|row_idx| {
            debug_assert_ne!(self.left_row_bounds[row_idx], usize::MAX);
            debug_assert_ne!(self.right_row_bounds[row_idx], 0);
        });
    }

    pub fn fill_rectangle(
        &mut self,
        target_start: usize,
        profile_start: usize,
        target_end: usize,
        profile_end: usize,
    ) {
        self.reuse(target_end);

        self.seq_start = target_start;
        self.seq_end = target_end;

        for row_idx in self.seq_start..=self.seq_end {
            self.left_row_bounds[row_idx] = profile_start;
            self.right_row_bounds[row_idx] = profile_end;
        }
    }

    pub fn valid(&self) -> bool {
        let mut prev_row_range = (
            self.left_row_bounds[self.seq_start],
            self.right_row_bounds[self.seq_start],
        );

        for row_idx in self.seq_start..=self.seq_end {
            let row_range = (
                self.left_row_bounds[row_idx],
                self.right_row_bounds[row_idx],
            );

            if row_range.0 > row_range.1 {
                return false;
            }

            // check if this range overlaps with the previous one
            if prev_row_range.0 > row_range.1 || prev_row_range.1 < row_range.0 {
                return false;
            }

            prev_row_range = row_range;
        }
        true
    }

    pub fn count_cells(&self) -> usize {
        let mut cells = 0;
        for target_idx in self.seq_start..=self.seq_end {
            cells += self.right_row_bounds[target_idx] - self.left_row_bounds[target_idx] + 1;
        }
        cells
    }

    pub fn dump(&self, out: &mut impl Write) -> Result<()> {
        for row_idx in self.seq_start..=self.seq_end {
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
        writeln!(f, "rows: {}-{}", self.seq_start, self.seq_end)?;
        for row_idx in self.seq_start..=self.seq_end {
            writeln!(
                f,
                "{}: {}-{}",
                row_idx, self.left_row_bounds[row_idx], self.right_row_bounds[row_idx]
            )?;
        }
        Ok(())
    }
}
