use anyhow::bail;

use crate::util::PrintMe;
use std::fmt::{Debug, Formatter};
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

pub struct Interval {
    pub start: usize,
    pub end: usize,
}

pub struct BoundingBox {
    target_start: usize,
    target_end: usize,
    profile_start: usize,
    profile_end: usize,
}

impl Default for BoundingBox {
    fn default() -> Self {
        Self {
            target_start: usize::MAX,
            target_end: 0,
            profile_start: usize::MAX,
            profile_end: 0,
        }
    }
}

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
    #[allow(dead_code)]
    fn new(
        left_target_idx: usize,
        left_profile_idx: usize,
        right_target_idx: usize,
        right_profile_idx: usize,
    ) -> Self {
        Self {
            left_target_idx,
            left_profile_idx,
            right_target_idx,
            right_profile_idx,
        }
    }

    pub fn reset(&mut self) {
        self.left_target_idx = 0;
        self.left_profile_idx = 1;
        self.right_target_idx = 1;
        self.right_profile_idx = 0;
    }

    pub fn is_default(&self) -> bool {
        self.left_target_idx == 0
            && self.left_profile_idx == 1
            && self.right_target_idx == 1
            && self.right_profile_idx == 0
    }

    pub fn is_single_cell(&self) -> bool {
        self.left_target_idx == self.right_target_idx
            && self.left_profile_idx == self.right_profile_idx
    }

    pub fn was_pruned(&self) -> bool {
        self.left_target_idx < self.right_target_idx
    }

    pub fn valid(&self) -> bool {
        // the computed anti-diagonal indices should match
        self.left_target_idx + self.left_profile_idx
            == self.right_target_idx + self.right_profile_idx
            // the left cell should be lower than the right cell
            && self.left_target_idx >= self.right_target_idx
            // the left cell should be left of the right cell
            && self.left_profile_idx <= self.right_profile_idx
    }

    pub fn idx(&self) -> usize {
        self.left_target_idx + self.left_profile_idx
    }

    pub fn cell_zip(&self) -> Zip<Rev<RangeInclusive<usize>>, RangeInclusive<usize>> {
        let target_range = (self.right_target_idx..=self.left_target_idx).rev();
        let profile_range = self.left_profile_idx..=self.right_profile_idx;
        target_range.zip(profile_range)
    }

    pub fn overlaps(&self, other: &Self) -> bool {
        if self.idx() != other.idx() {
            false
        } else {
            self.left_profile_idx <= other.right_profile_idx
                && self.right_profile_idx >= other.left_profile_idx
        }
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct AntiDiagonalBounds {
    pub bounds: Vec<AntiDiagonal>,
    pub target_length: usize,
    pub profile_length: usize,
    pub size: usize,
    pub min_anti_diagonal_idx: usize,
    pub max_anti_diagonal_idx: usize,
}

impl Default for AntiDiagonalBounds {
    fn default() -> Self {
        Self {
            bounds: vec![AntiDiagonal::default()],
            target_length: 0,
            profile_length: 0,
            size: 1,
            min_anti_diagonal_idx: 0,
            max_anti_diagonal_idx: 0,
        }
    }
}

impl AntiDiagonalBounds {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        let size = target_length + profile_length + 1;
        Self {
            bounds: vec![AntiDiagonal::default(); size],
            target_length,
            profile_length,
            size,
            min_anti_diagonal_idx: 0,
            max_anti_diagonal_idx: 0,
        }
    }

    pub fn resize(&mut self, new_size: usize) {
        if new_size > self.size {
            self.bounds.resize(new_size, AntiDiagonal::default());
            self.size = new_size;
        }
    }

    pub fn reset(&mut self) {
        self.bounds_mut().iter_mut().for_each(|bound| bound.reset());

        self.bounds
            .iter()
            .for_each(|b| debug_assert!(b.is_default()));
    }

    pub fn reuse(&mut self, target_length: usize, profile_length: usize) {
        self.reset();

        let new_size = target_length + profile_length + 1;
        self.resize(new_size);

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
        // **NOTE: these asserts always fail currently because
        //         of the access pattern during cloud search
        //
        // make sure the two cells are on the same anti-diagonal
        // debug_assert_eq!(
        //     left_target_idx + left_profile_idx,
        //     right_target_idx + right_profile_idx
        // );

        // make sure we are setting the anti-diagonal that we think we are
        //debug_assert_eq!(anti_diagonal_idx, left_target_idx + left_profile_idx);

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

    pub fn first(&self) -> &AntiDiagonal {
        self.get(self.min_anti_diagonal_idx)
    }

    pub fn last(&self) -> &AntiDiagonal {
        self.get(self.max_anti_diagonal_idx)
    }

    pub fn valid(&self) -> bool {
        self.bounds().iter().all(|bound| bound.valid())
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
    pub fn trim_wings(&mut self) -> anyhow::Result<()> {
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

            // the number of cells the right bound
            // is above the previous right bound
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

            // the number of cells the left bound
            // is below the previous left bound
            let left_distance =
                (previous_bound.left_profile_idx).saturating_sub(current_bound.left_profile_idx);

            self.set(
                anti_diagonal_idx,
                current_bound.left_target_idx - left_distance,
                current_bound.left_profile_idx + left_distance,
                current_bound.right_target_idx + right_distance,
                current_bound.right_profile_idx - right_distance,
            );

            if !self.get(anti_diagonal_idx).valid() {
                bail!("cloud trimming failed");
            }
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

            // the number of cells the right bound
            // is above the previous right bound
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

            // the number of cells the left bound
            // is below the previous left bound
            let left_distance =
                (current_bound.left_target_idx).saturating_sub(previous_bound.left_target_idx);

            self.set(
                anti_diagonal_idx,
                current_bound.left_target_idx - left_distance,
                current_bound.left_profile_idx + left_distance,
                current_bound.right_target_idx + right_distance,
                current_bound.right_profile_idx - right_distance,
            );

            if !self.get(anti_diagonal_idx).valid() {
                bail!("cloud trimming failed");
            }
        }

        Ok(())
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
            let first_bound = self.first();
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
            let last_bound = self.last();
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

        self.min_anti_diagonal_idx = first_anti_diagonal_idx - left_distance;
        self.max_anti_diagonal_idx = last_anti_diagonal_idx + right_distance;
    }

    pub fn advance_forward(&mut self) {
        let last_bound = self.last();

        let next_idx = self.max_anti_diagonal_idx + 1;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move past the target length, and 0 otherwise
        let profile_addend = ((last_bound.left_target_idx + 1) > self.target_length) as usize;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move past the profile length, and 0 otherwise
        let target_addend = ((last_bound.right_profile_idx + 1) > self.profile_length) as usize;

        self.set(
            next_idx,
            // prevent moving past the target length
            (last_bound.left_target_idx + 1).min(self.target_length),
            // if we hit the target length, then the
            // addend is 1 and we'll move right a cell
            last_bound.left_profile_idx + profile_addend,
            // if we hit the profile length, then the
            // addend is 1 and we'll drop down a cell
            last_bound.right_target_idx + target_addend,
            // prevent moving past the profile length
            (last_bound.right_profile_idx + 1).min(self.profile_length),
        );

        self.max_anti_diagonal_idx = next_idx;
    }

    pub fn advance_reverse(&mut self) {
        let first_bound = self.first();

        let next_idx = self.min_anti_diagonal_idx - 1;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move before target index 1, and 0 otherwise
        let profile_subtrahend = ((first_bound.right_target_idx.saturating_sub(1)) < 1) as usize;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move before profle index 1, and 0 otherwise
        let target_subtrahend = ((first_bound.left_profile_idx.saturating_sub(1)) < 1) as usize;

        self.set(
            next_idx,
            // if we hit the profile index 1, then the
            // subtrahend is 1 and we'll move up a cell
            first_bound.left_target_idx - target_subtrahend,
            // prevent moving before profile index 1
            (first_bound.left_profile_idx - 1).max(1),
            // prevent moving before target index 1
            (first_bound.right_target_idx - 1).max(1),
            // if we hit the target index 1, then the
            // subtrahend is 1 and we'll move left a cell
            first_bound.right_profile_idx - profile_subtrahend,
        );

        self.min_anti_diagonal_idx = next_idx;
    }

    pub fn fill_rectangle(
        &mut self,
        target_start: usize,
        profile_start: usize,
        target_end: usize,
        profile_end: usize,
    ) {
        self.reset();

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

    pub fn bounding_box(&self) -> BoundingBox {
        let mut bbox = BoundingBox::default();

        self.bounds().iter().for_each(|bound| {
            bbox.target_start = bbox.target_start.min(bound.left_target_idx);
            bbox.target_end = bbox.target_start.max(bound.right_target_idx);
            bbox.profile_start = bbox.profile_start.min(bound.left_profile_idx);
            bbox.profile_end = bbox.profile_start.max(bound.right_profile_idx);
        });

        bbox
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

    pub fn intersection(b1: &Self, b2: &Self) -> Option<Interval> {
        let overlap_start = b1.min_anti_diagonal_idx.max(b2.min_anti_diagonal_idx);
        let overlap_end = b1.max_anti_diagonal_idx.min(b2.max_anti_diagonal_idx);

        let anti_diagonal_ranges_overlap = overlap_start <= overlap_end;

        if !anti_diagonal_ranges_overlap {
            return None;
        }

        let mut start = None;
        for anti_diagonal_idx in overlap_start..=overlap_end {
            let a_bound = b1.get(anti_diagonal_idx);
            let b_bound = b2.get(anti_diagonal_idx);

            if a_bound.overlaps(b_bound) {
                start = Some(anti_diagonal_idx);
                break;
            }
        }

        let mut end = None;
        for anti_diagonal_idx in (overlap_start..=overlap_end).rev() {
            let a_bound = b1.get(anti_diagonal_idx);
            let b_bound = b2.get(anti_diagonal_idx);

            if a_bound.overlaps(b_bound) {
                end = Some(anti_diagonal_idx);
                break;
            }
        }

        match (start, end) {
            (Some(start), Some(end)) => Some(Interval { start, end }),
            (None, None) => None,
            (_, _) => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_advance_forward() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(2, 1, 1, 1, 1);
        bounds.min_anti_diagonal_idx = 2;
        bounds.max_anti_diagonal_idx = 2;

        (0..18).for_each(|_| bounds.advance_forward());

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        target_bounds.fill_rectangle(1, 1, 10, 10);

        assert!(bounds == target_bounds);
    }

    #[test]
    pub fn test_advance_reverse() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(20, 10, 10, 10, 10);
        bounds.min_anti_diagonal_idx = 20;
        bounds.max_anti_diagonal_idx = 20;

        (0..18).for_each(|_| bounds.advance_reverse());

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        target_bounds.fill_rectangle(1, 1, 10, 10);

        assert!(bounds == target_bounds);
    }

    #[test]
    pub fn test_fill_square() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.fill_rectangle(3, 3, 6, 6);

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        target_bounds.set(6, 3, 3, 3, 3);
        target_bounds.set(7, 4, 3, 3, 4);
        target_bounds.set(8, 5, 3, 3, 5);
        target_bounds.set(9, 6, 3, 3, 6);
        target_bounds.set(10, 6, 4, 4, 6);
        target_bounds.set(11, 6, 5, 5, 6);
        target_bounds.set(12, 6, 6, 6, 6);
        target_bounds.min_anti_diagonal_idx = 6;
        target_bounds.max_anti_diagonal_idx = 12;

        assert!(bounds == target_bounds);
    }

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

        bounds.min_anti_diagonal_idx = 4;
        bounds.max_anti_diagonal_idx = 16;

        let mut target_bounds = bounds.clone();
        target_bounds.set(8, 5, 3, 3, 5);
        target_bounds.set(12, 7, 5, 5, 7);

        bounds.trim_wings();
        assert!(bounds == target_bounds);

        // re-run it to make sure it does
        // nothing to already-trimmed bounds
        bounds.trim_wings();
        assert!(bounds == target_bounds);
    }

    #[test]
    fn test_square_corners() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(9, 6, 3, 3, 6);
        bounds.set(10, 6, 4, 4, 6);
        bounds.set(11, 7, 4, 4, 7);

        bounds.min_anti_diagonal_idx = 9;
        bounds.max_anti_diagonal_idx = 11;

        let mut target_bounds = bounds.clone();
        target_bounds.set(6, 3, 3, 3, 3);
        target_bounds.set(7, 4, 3, 3, 4);
        target_bounds.set(8, 5, 3, 3, 5);

        target_bounds.set(12, 7, 5, 5, 7);
        target_bounds.set(13, 7, 6, 6, 7);
        target_bounds.set(14, 7, 7, 7, 7);

        target_bounds.min_anti_diagonal_idx = 6;
        target_bounds.max_anti_diagonal_idx = 14;

        bounds.square_corners();
        assert!(bounds == target_bounds);

        // re-run it to make sure it does
        // nothing to already-squared bounds
        bounds.square_corners();
        assert!(bounds == target_bounds);
    }

    #[test]
    fn test_merge_join() {
        let mut forward_bounds = AntiDiagonalBounds::new(10, 10);
        forward_bounds.set(8, 4, 4, 4, 4);
        forward_bounds.set(9, 5, 4, 4, 5);
        forward_bounds.set(10, 6, 4, 4, 6);
        forward_bounds.set(11, 7, 4, 4, 7);
        forward_bounds.set(12, 8, 4, 4, 8);
        forward_bounds.set(13, 8, 5, 5, 8);
        forward_bounds.set(14, 8, 6, 6, 8);

        forward_bounds.min_anti_diagonal_idx = 8;
        forward_bounds.max_anti_diagonal_idx = 14;

        let mut reverse_bounds = AntiDiagonalBounds::new(10, 10);
        reverse_bounds.set(6, 4, 2, 2, 4);
        reverse_bounds.set(7, 5, 2, 2, 5);
        reverse_bounds.set(8, 6, 2, 2, 6);
        reverse_bounds.set(9, 6, 3, 3, 6);
        reverse_bounds.set(10, 6, 4, 4, 6);
        reverse_bounds.set(11, 6, 5, 5, 6);
        reverse_bounds.set(12, 6, 6, 6, 6);

        reverse_bounds.min_anti_diagonal_idx = 6;
        reverse_bounds.max_anti_diagonal_idx = 12;

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        target_bounds.set(6, 4, 2, 2, 4);
        target_bounds.set(7, 5, 2, 2, 5);
        target_bounds.set(8, 6, 2, 2, 6);
        target_bounds.set(9, 6, 3, 3, 6);
        target_bounds.set(10, 6, 4, 4, 6);
        target_bounds.set(11, 7, 4, 4, 7);
        target_bounds.set(12, 8, 4, 4, 8);
        target_bounds.set(13, 8, 5, 5, 8);
        target_bounds.set(14, 8, 6, 6, 8);

        target_bounds.min_anti_diagonal_idx = 6;
        target_bounds.max_anti_diagonal_idx = 14;

        AntiDiagonalBounds::join_merge(&mut forward_bounds, &reverse_bounds);
        assert!(forward_bounds == target_bounds);

        // re-run it to make sure that
        // re-joining the bounds does nothing
        AntiDiagonalBounds::join_merge(&mut forward_bounds, &reverse_bounds);
        assert!(forward_bounds == target_bounds);
    }

    #[test]
    fn test_anti_diagonal_overlaps() {
        // different anti-diagonals
        let a = AntiDiagonal::new(9, 1, 1, 9);
        let b = AntiDiagonal::new(9, 2, 2, 9);
        assert!(!a.overlaps(&b));

        // same antidiagonal ovlerap
        let a = AntiDiagonal::new(9, 1, 1, 9);
        let b = AntiDiagonal::new(9, 1, 1, 9);
        assert!(a.overlaps(&b));

        // single cell overlap
        let a = AntiDiagonal::new(9, 1, 5, 5);
        let b = AntiDiagonal::new(5, 5, 1, 9);
        assert!(a.overlaps(&b));

        // single cell difference
        let a = AntiDiagonal::new(9, 1, 6, 4);
        let b = AntiDiagonal::new(5, 5, 1, 9);
        assert!(!a.overlaps(&b));

        // several cell overlap
        let a = AntiDiagonal::new(9, 1, 3, 7);
        let b = AntiDiagonal::new(7, 3, 1, 9);
        assert!(a.overlaps(&b));

        // several cell difference
        let a = AntiDiagonal::new(9, 1, 7, 3);
        let b = AntiDiagonal::new(3, 7, 1, 9);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_anti_diagonal_bounds_intersects() {
        // no overlap, not touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 4, 4);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(6, 6, 10, 10);

        assert!(a.intersects(&a));
        assert!(b.intersects(&b));
        assert!(!a.intersects(&b));
        assert!(!b.intersects(&a));

        // no overlap, touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 5, 5);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(6, 6, 10, 10);

        assert!(a.intersects(&a));
        assert!(b.intersects(&b));
        assert!(!a.intersects(&b));
        assert!(!b.intersects(&a));

        // 1 cell overlap
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 5, 5);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(5, 5, 10, 10);

        assert!(a.intersects(&a));
        assert!(b.intersects(&b));
        assert!(a.intersects(&b));
        assert!(b.intersects(&a));

        // many cell overlap
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 8, 8);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(2, 2, 10, 10);

        assert!(a.intersects(&a));
        assert!(b.intersects(&b));
        assert!(a.intersects(&b));
        assert!(b.intersects(&a));

        // anti-diagonal overlap, not touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 8, 4);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(1, 5, 8, 9);

        assert!(a.intersects(&a));
        assert!(b.intersects(&b));
        assert!(!a.intersects(&b));
        assert!(!b.intersects(&a));
    }
}
