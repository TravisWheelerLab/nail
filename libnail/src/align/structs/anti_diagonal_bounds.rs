use crate::util::PrintMe;
use std::fmt::{Debug, Formatter};
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

#[derive(Clone, PartialEq)]
pub struct AntiDiagonal {
    pub left_target_idx: usize,
    pub left_profile_idx: usize,
    pub right_target_idx: usize,
    pub right_profile_idx: usize,
}

impl Default for AntiDiagonal {
    fn default() -> Self {
        // set the default to a simple invalid bound
        // (the left bound is on the right, and vice versa)
        AntiDiagonal {
            left_target_idx: 0,
            left_profile_idx: 1,
            right_target_idx: 1,
            right_profile_idx: 0,
        }
    }
}

impl Debug for AntiDiagonal {
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

impl PrintMe for AntiDiagonal {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl AntiDiagonal {
    pub fn was_pruned(&self) -> bool {
        self.left_target_idx < self.right_target_idx
    }

    pub fn idx(&self) -> usize {
        self.left_target_idx + self.left_profile_idx
    }

    pub fn cell_zip(&self) -> Zip<Rev<RangeInclusive<usize>>, RangeInclusive<usize>> {
        let target_range = (self.right_target_idx..=self.left_target_idx).rev();
        let profile_range = self.left_profile_idx..=self.right_profile_idx;
        target_range.zip(profile_range)
    }
}

#[derive(Default, Clone)]
pub struct AntiDiagonalBounds {
    pub bounds: Vec<AntiDiagonal>,
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
            bounds: vec![AntiDiagonal::default(); size],
            target_length,
            profile_length,
            size,
            min_anti_diagonal_idx: size,
            max_anti_diagonal_idx: 0,
        }
    }

    pub fn resize(&mut self, new_size: usize) {
        self.bounds.resize(new_size, AntiDiagonal::default());
        self.size = new_size;
    }

    pub fn reuse(&mut self, target_length: usize, profile_length: usize) {
        let new_size = target_length + profile_length + 1;
        if new_size > self.size {
            self.resize(new_size);
        }

        for bound in self.bounds_mut() {
            bound.left_target_idx = 0;
            bound.left_profile_idx = 1;
            bound.right_target_idx = 1;
            bound.right_profile_idx = 0;
        }

        self.min_anti_diagonal_idx = new_size;
        self.max_anti_diagonal_idx = 0;

        self.target_length = target_length;
        self.profile_length = profile_length;
    }

    pub fn set(
        &mut self,
        anti_diagonal_idx: usize,
        left_target_idx: usize,
        left_profile_idx: usize,
        right_target_idx: usize,
        right_profile_idx: usize,
    ) {
        // make sure the two cells are on the same anti-diagonal
        debug_assert_eq!(
            left_target_idx + left_profile_idx,
            right_target_idx + right_profile_idx
        );

        // make sure we are setting the anti-diagonal that we think we are
        debug_assert_eq!(anti_diagonal_idx, left_target_idx + left_profile_idx);

        self.min_anti_diagonal_idx = self.min_anti_diagonal_idx.min(anti_diagonal_idx);
        self.max_anti_diagonal_idx = self.max_anti_diagonal_idx.max(anti_diagonal_idx);

        let bound = &mut self.bounds[anti_diagonal_idx];
        bound.left_target_idx = left_target_idx;
        bound.left_profile_idx = left_profile_idx;
        bound.right_target_idx = right_target_idx;
        bound.right_profile_idx = right_profile_idx;
    }

    pub fn get(&self, idx: usize) -> &AntiDiagonal {
        &self.bounds[idx]
    }

    pub fn get_mut(&mut self, idx: usize) -> &mut AntiDiagonal {
        &mut self.bounds[idx]
    }

    pub fn get_first(&self) -> &AntiDiagonal {
        self.get(self.min_anti_diagonal_idx)
    }

    pub fn get_last(&self) -> &AntiDiagonal {
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

    pub fn bounds(&self) -> &[AntiDiagonal] {
        &self.bounds[self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx]
    }

    pub fn bounds_mut(&mut self) -> &mut [AntiDiagonal] {
        &mut self.bounds[self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx]
    }

    /// This removes all of the protruding regions in the cloud that are
    /// unreachable from a traceback that traverses the entire cloud.
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
    pub fn square_corners(&mut self) {
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
                first_bound.idx(),
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
                last_bound.idx(),
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

            debug_assert_eq!(forward_bound.idx(), backward_bound.idx());

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_trim_wings() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(4, 2, 2, 2, 2);
        bounds.set(5, 3, 2, 2, 3);
        bounds.set(6, 4, 2, 2, 4);
        bounds.set(7, 5, 2, 2, 5);
        bounds.set(8, 6, 2, 2, 6);
        bounds.set(9, 5, 4, 4, 5);
        bounds.set(10, 6, 4, 4, 6);
        bounds.set(11, 6, 5, 5, 6);
        bounds.set(12, 8, 4, 4, 8);
        bounds.set(13, 8, 5, 5, 8);
        bounds.set(14, 8, 6, 6, 8);
        bounds.set(15, 8, 7, 7, 8);
        bounds.set(16, 8, 8, 8, 8);

        let mut trimmed_bounds = bounds.clone();

        trimmed_bounds.trim_wings();

        assert!(trimmed_bounds.valid());

        let invalid_range = (0..=3).chain(17..=20);
        invalid_range
            .map(|idx| trimmed_bounds.get(idx))
            .for_each(|b| assert_eq!(*b, AntiDiagonal::default()));

        let unchanged_range = (4..=7).chain(9..=11).chain(13..=16);
        unchanged_range
            .map(|idx| (trimmed_bounds.get(idx), bounds.get(idx)))
            .for_each(|(b1, b2)| assert_eq!(*b1, *b2));

        let b = trimmed_bounds.get(8);
        assert_eq!(b.left_target_idx, 5);
        assert_eq!(b.left_profile_idx, 3);
        assert_eq!(b.right_target_idx, 3);
        assert_eq!(b.right_profile_idx, 5);

        let b = trimmed_bounds.get(12);
        assert_eq!(b.left_target_idx, 7);
        assert_eq!(b.left_profile_idx, 5);
        assert_eq!(b.right_target_idx, 5);
        assert_eq!(b.right_profile_idx, 7);

        let mut retrimmed_bounds = trimmed_bounds.clone();
        retrimmed_bounds.trim_wings();

        trimmed_bounds
            .bounds
            .into_iter()
            .zip(retrimmed_bounds.bounds)
            .for_each(|(b1, b2)| {
                assert_eq!(b1, b2);
            })
    }

    #[test]
    fn test_square_corners() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(9, 6, 3, 3, 6);
        bounds.set(10, 6, 4, 4, 6);
        bounds.set(11, 7, 4, 4, 7);

        let mut squared_bounds = bounds.clone();
        squared_bounds.square_corners();

        let invalid_range = (0..=5).chain(15..=20);
        invalid_range
            .map(|idx| squared_bounds.get(idx))
            .for_each(|b| assert_eq!(*b, AntiDiagonal::default()));

        let unchanged_range = 9..=11;
        unchanged_range
            .map(|idx| (squared_bounds.get(idx), bounds.get(idx)))
            .for_each(|(b1, b2)| assert_eq!(b1, b2));

        let b = squared_bounds.get(6);
        assert_eq!(b.left_target_idx, 3);
        assert_eq!(b.left_profile_idx, 3);
        assert_eq!(b.right_target_idx, 3);
        assert_eq!(b.right_profile_idx, 3);

        let b = squared_bounds.get(7);
        assert_eq!(b.left_target_idx, 4);
        assert_eq!(b.left_profile_idx, 3);
        assert_eq!(b.right_target_idx, 3);
        assert_eq!(b.right_profile_idx, 4);

        let b = squared_bounds.get(8);
        assert_eq!(b.left_target_idx, 5);
        assert_eq!(b.left_profile_idx, 3);
        assert_eq!(b.right_target_idx, 3);
        assert_eq!(b.right_profile_idx, 5);

        let b = squared_bounds.get(12);
        assert_eq!(b.left_target_idx, 7);
        assert_eq!(b.left_profile_idx, 5);
        assert_eq!(b.right_target_idx, 5);
        assert_eq!(b.right_profile_idx, 7);

        let b = squared_bounds.get(13);
        assert_eq!(b.left_target_idx, 7);
        assert_eq!(b.left_profile_idx, 6);
        assert_eq!(b.right_target_idx, 6);
        assert_eq!(b.right_profile_idx, 7);

        let b = squared_bounds.get(14);
        assert_eq!(b.left_target_idx, 7);
        assert_eq!(b.left_profile_idx, 7);
        assert_eq!(b.right_target_idx, 7);
        assert_eq!(b.right_profile_idx, 7);
    }
}
