use crate::util::PrintMe;
use std::fmt::{Display, Formatter};
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

impl Display for CloudBound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{},{} : {},{}",
            self.left_target_idx,
            self.left_profile_idx,
            self.right_target_idx,
            self.right_profile_idx
        )
    }
}

impl PrintMe for CloudBound {
    fn print(&self) {
        println!("{}", self);
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

#[derive(Default, Clone)]
pub struct AntiDiagonalBounds {
    pub bounds: Vec<CloudBound>,
    pub target_length: usize,
    pub profile_length: usize,
    pub size: usize,
    pub min_anti_diagonal_idx: usize,
    pub max_anti_diagonal_idx: usize,
}

impl AntiDiagonalBounds {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        let size = target_length + profile_length + 1;
        Self {
            bounds: vec![CloudBound::default(); size],
            target_length,
            profile_length,
            size,
            min_anti_diagonal_idx: size,
            max_anti_diagonal_idx: 0,
        }
    }

    pub fn resize(&mut self, new_size: usize) {
        self.bounds.resize(new_size, CloudBound::default());
        self.size = new_size;
    }

    pub fn reuse(&mut self, target_length: usize, profile_length: usize) {
        let new_size = target_length + profile_length + 1;
        if new_size > self.size {
            self.resize(new_size);
        }

        self.min_anti_diagonal_idx = new_size;
        self.max_anti_diagonal_idx = 0;

        self.target_length = target_length;
        self.profile_length = profile_length;

        // TODO: think about this
        for bound in self.bounds.iter_mut() {
            bound.left_target_idx = 0;
            bound.left_profile_idx = 1;
            bound.right_target_idx = 1;
            bound.right_profile_idx = 0;
        }
    }

    pub fn fill(&mut self) {
        todo!()
    }

    pub fn set(
        &mut self,
        // TODO: remove this
        anti_diagonal_idx: usize,
        left_target_idx: usize,
        left_profile_idx: usize,
        right_target_idx: usize,
        right_profile_idx: usize,
    ) {
        // TODO: I don't think there's much point to setting these here
        self.min_anti_diagonal_idx = self.min_anti_diagonal_idx.min(anti_diagonal_idx);
        self.max_anti_diagonal_idx = self.max_anti_diagonal_idx.max(anti_diagonal_idx);

        debug_assert_eq!(
            left_target_idx + left_profile_idx,
            right_target_idx + right_profile_idx
        );

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

    pub fn get_first(&self) -> &CloudBound {
        self.get(self.min_anti_diagonal_idx)
    }

    pub fn get_last(&self) -> &CloudBound {
        self.get(self.max_anti_diagonal_idx)
    }

    pub fn valid(&self) -> bool {
        for bound in self.bounds() {
            if bound.was_pruned() {
                return false;
            }
        }
        true
    }

    /// Get the total number of cells that exist within the cloud boundaries.
    pub fn cloud_size(&self) -> usize {
        let mut cloud_size = 0usize;
        for bound in self.bounds[self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx].iter() {
            cloud_size += bound.left_target_idx - bound.right_target_idx;
        }
        cloud_size
    }

    /// Get the number of anti-diagonals defined in the cloud.
    pub fn num_anti_diagonals(&self) -> usize {
        self.max_anti_diagonal_idx - self.min_anti_diagonal_idx + 1
    }

    /// This removes all of the protruding regions in the cloud that are unreachable from a traceback
    /// that traverses the entire cloud.
    pub fn trim_wings(&mut self) {
        for anti_diagonal_idx in self.min_anti_diagonal_idx + 1..=self.max_anti_diagonal_idx {
            let previous_bound = self.get(anti_diagonal_idx - 1);
            let current_bound = self.get(anti_diagonal_idx);

            // right wings identifiable from a forward pass look like:
            //            c
            //         c  c
            //   c  c  c  c
            //   c  c  c  c
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //
            let right_distance = previous_bound
                .right_target_idx
                .saturating_sub(current_bound.right_target_idx);

            // left wings identifiable from a forward pass look like:
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //         c  c  c  c
            //      c  c  c  c  c
            //   c  c  c  c  c  c
            //
            let left_distance =
                (previous_bound.left_profile_idx).saturating_sub(current_bound.left_profile_idx);

            self.set(
                anti_diagonal_idx,
                current_bound.left_target_idx - left_distance,
                current_bound.left_profile_idx + left_distance,
                current_bound.right_target_idx + right_distance,
                current_bound.right_profile_idx - right_distance,
            )
        }

        for anti_diagonal_idx in (self.min_anti_diagonal_idx..self.max_anti_diagonal_idx).rev() {
            let previous_bound = self.get(anti_diagonal_idx + 1);
            let current_bound = self.get(anti_diagonal_idx);

            // right wings identifiable from a backward pass look like:
            //   c  c  c  c  c  c
            //   c  c  c  c  c
            //   c  c  c  c
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //
            let right_distance = current_bound
                .right_profile_idx
                .saturating_sub(previous_bound.right_profile_idx);

            // left wings identifiable from a backward pass look like:
            //   c  c  c  c  c  c
            //   c  c  c  c  c  c
            //         c  c  c  c
            //         c  c  c  c
            //         c  c
            //         c
            //
            let left_distance =
                (current_bound.left_target_idx).saturating_sub(previous_bound.left_target_idx);

            self.set(
                anti_diagonal_idx,
                current_bound.left_target_idx - left_distance,
                current_bound.left_profile_idx + left_distance,
                current_bound.right_target_idx + right_distance,
                current_bound.right_profile_idx - right_distance,
            )
        }
    }

    /// Add a square block of bounds to the group with side length of `size`.
    ///
    /// This is currently intended for debugging.
    pub fn bound_block(&mut self, target_start: usize, profile_start: usize, size: usize) {
        let idx = target_start + profile_start;
        for i in 0..size {
            self.set(
                idx + i,
                target_start + i,
                profile_start,
                target_start,
                profile_start + i,
            );
        }
        let idx = target_start + profile_start + size;
        for i in 0..=size {
            self.set(
                idx + i,
                target_start + size,
                profile_start + i,
                target_start + i,
                profile_start + size,
            );
        }
    }

    pub fn bounds(&self) -> &[CloudBound] {
        &self.bounds[self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx]
    }

    pub fn square_corners(&mut self) {
        // TODO: refactor this nastiness once I've
        //       refactored the AD bounds struct
        // this is a bit of jumping through hoops
        // to appease the borrow checker
        let (
            left_distance,
            first_anti_diagonal_idx,
            first_left_target_idx,
            first_left_profile_idx,
            first_right_target_idx,
            first_right_profile_idx,
        ) = {
            let first_bound = self.get_first();
            (
                first_bound.left_target_idx - first_bound.right_target_idx,
                first_bound.anti_diagonal_idx(),
                first_bound.left_target_idx,
                first_bound.left_profile_idx,
                first_bound.right_target_idx,
                first_bound.right_profile_idx,
            )
        };

        (1..=left_distance)
            .map(|offset| (offset, first_anti_diagonal_idx - offset))
            .for_each(|(offset, anti_diagonal_idx)| {
                self.set(
                    anti_diagonal_idx,
                    first_left_target_idx - offset,
                    first_left_profile_idx,
                    first_right_target_idx,
                    first_right_profile_idx - offset,
                );
            });
        let (
            right_distance,
            last_anti_diagonal_idx,
            last_left_target_idx,
            last_left_profile_idx,
            last_right_target_idx,
            last_right_profile_idx,
        ) = {
            let last_bound = self.get_last();
            (
                last_bound.left_target_idx - last_bound.right_target_idx,
                last_bound.anti_diagonal_idx(),
                last_bound.left_target_idx,
                last_bound.left_profile_idx,
                last_bound.right_target_idx,
                last_bound.right_profile_idx,
            )
        };

        (1..=right_distance)
            .map(|offset| (offset, last_anti_diagonal_idx + offset))
            .for_each(|(offset, anti_diagonal_idx)| {
                self.set(
                    anti_diagonal_idx,
                    last_left_target_idx,
                    last_left_profile_idx + offset,
                    last_right_target_idx + offset,
                    last_right_profile_idx,
                );
            });
    }

    pub fn fill_rectangle(
        &mut self,
        target_start: usize,
        profile_start: usize,
        target_end: usize,
        profile_end: usize,
    ) {
        let target_distance = target_end - target_start;
        let profile_distance = profile_end - profile_start;

        let anti_diagonal_start = target_start + profile_start;
        let anti_diagonal_end = target_end + profile_end;

        for idx in anti_diagonal_start..=anti_diagonal_end {
            // relative_idx is the number of antidiagonals we have moved
            // forward relative to the starting antidiagonal index
            let relative_idx = idx - anti_diagonal_start;
            let bound = self.get_mut(idx);

            bound.left_target_idx = target_end.min(target_start + relative_idx);
            bound.left_profile_idx = profile_start + relative_idx.saturating_sub(target_distance);

            bound.right_target_idx = target_start + relative_idx.saturating_sub(profile_distance);
            bound.right_profile_idx = profile_end.min(profile_start + relative_idx);

            debug_assert_eq!(
                bound.left_target_idx + bound.left_profile_idx,
                bound.right_target_idx + bound.right_profile_idx
            );
        }

        self.min_anti_diagonal_idx = anti_diagonal_start;
        self.max_anti_diagonal_idx = anti_diagonal_end;
    }

    pub fn join_merge(
        forward_bounds: &mut AntiDiagonalBounds,
        backward_bounds: &AntiDiagonalBounds,
    ) {
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

            debug_assert_eq!(
                forward_bound.anti_diagonal_idx(),
                backward_bound.anti_diagonal_idx()
            );

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
    }
}
