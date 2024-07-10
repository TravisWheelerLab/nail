use anyhow::bail;

use crate::util::PrintMe;
use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub fn len(&self) -> usize {
        self.end - self.start + 1
    }
}

#[derive(Debug, PartialEq)]
pub enum Relationship {
    Disjoint(Interval),
    Intersecting(Interval),
}

impl Relationship {
    pub fn interval(&self) -> Interval {
        match self {
            Relationship::Disjoint(interval) | Relationship::Intersecting(interval) => *interval,
        }
    }
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

    pub fn intersects(&self, other: &Self) -> bool {
        if self.idx() != other.idx() {
            false
        } else {
            self.left_profile_idx <= other.right_profile_idx
                && self.right_profile_idx >= other.left_profile_idx
        }
    }

    pub fn merge(&mut self, other: &Self) {
        self.left_target_idx = self.left_target_idx.max(other.left_target_idx);
        self.left_profile_idx = self.left_profile_idx.min(other.left_profile_idx);
        self.right_target_idx = self.right_target_idx.min(other.right_target_idx);
        self.right_profile_idx = self.right_profile_idx.max(other.right_profile_idx);
    }

    pub fn replace_with(&mut self, other: &Self) {
        self.left_target_idx = other.left_target_idx;
        self.left_profile_idx = other.left_profile_idx;
        self.right_target_idx = other.right_target_idx;
        self.right_profile_idx = other.right_profile_idx;
    }

    pub fn grow_up(&mut self, target_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.right_target_idx = self.right_target_idx.min(target_idx).max(1);
        self.right_profile_idx = idx - self.right_target_idx;

        debug_assert_eq!(idx, self.idx())
    }

    pub fn grow_down(&mut self, target_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.left_target_idx = self.left_target_idx.max(target_idx);
        self.left_profile_idx = idx - self.left_target_idx;

        debug_assert_eq!(idx, self.idx())
    }

    pub fn grow_left(&mut self, profile_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.left_profile_idx = self.left_profile_idx.min(profile_idx).max(1);
        self.left_target_idx = idx - self.left_profile_idx;

        debug_assert_eq!(idx, self.idx())
    }

    pub fn grow_right(&mut self, profile_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.right_profile_idx = self.right_profile_idx.max(profile_idx);
        self.right_target_idx = idx - self.right_profile_idx;

        debug_assert_eq!(idx, self.idx())
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

    pub fn left_offset_to(&mut self, other: &Self) -> Ordering {
        self.min_anti_diagonal_idx.cmp(&other.min_anti_diagonal_idx)
    }

    pub fn right_offset_to(&mut self, other: &Self) -> Ordering {
        self.max_anti_diagonal_idx.cmp(&other.max_anti_diagonal_idx)
    }

    pub fn merge(&mut self, other: &Self) {
        let relationship = self.anti_diagonal_relationship(other);
        let interval = relationship.interval();

        // if 'self' is right-offset from 'other', we need
        // to pull in the remaining 'other' bounds
        if let Ordering::Greater = self.left_offset_to(other) {
            let self_slice = &mut self.bounds[other.min_anti_diagonal_idx..interval.start];
            let other_slice = &other.bounds[other.min_anti_diagonal_idx..interval.start];

            self_slice
                .iter_mut()
                .zip(other_slice)
                .for_each(|(self_bound, other_bound)| {
                    self_bound.replace_with(other_bound);
                });

            self.min_anti_diagonal_idx = other.min_anti_diagonal_idx;
        }

        // if 'self' is left-offset from 'other', we need
        // to pull in the remaining 'other' bounds
        if let Ordering::Less = self.right_offset_to(other) {
            let self_slice = &mut self.bounds[(interval.end + 1)..=other.max_anti_diagonal_idx];
            let other_slice = &other.bounds[(interval.end + 1)..=other.max_anti_diagonal_idx];

            self_slice
                .iter_mut()
                .zip(other_slice)
                .for_each(|(self_bound, other_bound)| {
                    self_bound.replace_with(other_bound);
                });

            self.max_anti_diagonal_idx = other.max_anti_diagonal_idx;
        }

        match relationship {
            Relationship::Intersecting(_) => {
                // merge along the overlapping anti-diagonals
                let self_slice = &mut self.bounds[interval.start..=interval.end];
                let other_slice = &other.bounds[interval.start..=interval.end];

                self_slice
                    .iter_mut()
                    .zip(other_slice)
                    .for_each(|(self_bound, other_bound)| self_bound.merge(other_bound));

                let min_target_idx = self.get(interval.start).right_target_idx;
                let min_profile_idx = self.get(interval.start).left_profile_idx;

                let self_slice = &mut self.bounds[self.min_anti_diagonal_idx..interval.start];

                self_slice.iter_mut().for_each(|bound| {
                    bound.grow_up(min_target_idx);
                    bound.grow_left(min_profile_idx);

                    debug_assert!(bound.right_target_idx >= 1);
                    debug_assert!(bound.left_profile_idx >= 1);
                });

                let max_target_idx = self.get(interval.end).left_target_idx;
                let max_profile_idx = self.get(interval.end).right_profile_idx;

                let self_slice = &mut self.bounds[(interval.end + 1)..=self.max_anti_diagonal_idx];

                self_slice.iter_mut().for_each(|bound| {
                    bound.grow_right(max_profile_idx);
                    bound.grow_down(max_target_idx);

                    debug_assert!(bound.left_target_idx <= self.target_length);
                    debug_assert!(bound.left_profile_idx <= self.profile_length);
                });
            }
            Relationship::Disjoint(_) => {
                let target_start = self.get(interval.start - 1).right_target_idx;
                let profile_start = self.get(interval.start - 1).left_profile_idx;
                let target_end = self.get(interval.end + 1).left_target_idx;
                let profile_end = self.get(interval.end + 1).right_profile_idx;

                let target_distance = target_end - target_start;
                let profile_distance = profile_end - profile_start;
                let left_corner_idx = target_start + profile_start;

                self.bounds
                    .iter_mut()
                    .enumerate()
                    .skip(interval.start)
                    .take(interval.len())
                    .for_each(|(idx, bound)| {
                        let relative_idx = idx - left_corner_idx;

                        bound.left_target_idx = target_end.min(target_start + relative_idx);
                        bound.left_profile_idx =
                            profile_start + relative_idx.saturating_sub(target_distance);
                        bound.right_target_idx =
                            target_start + relative_idx.saturating_sub(profile_distance);
                        bound.right_profile_idx = profile_end.min(profile_start + relative_idx);
                    });
            }
        }
    }

    pub fn anti_diagonal_relationship(&self, other: &Self) -> Relationship {
        let overlap_start = self.min_anti_diagonal_idx.max(other.min_anti_diagonal_idx);
        let overlap_end = self.max_anti_diagonal_idx.min(other.max_anti_diagonal_idx);

        let anti_diagonal_ranges_overlap = overlap_start <= overlap_end;

        if anti_diagonal_ranges_overlap {
            Relationship::Intersecting(Interval {
                start: overlap_start,
                end: overlap_end,
            })
        } else {
            Relationship::Disjoint(Interval {
                start: overlap_end + 1,
                end: overlap_start - 1,
            })
        }
    }

    pub fn cloud_relationship(&self, other: &Self) -> Relationship {
        match self.anti_diagonal_relationship(other) {
            // if the AD relationship is intersecting, we
            // still need to check if any of the ADs on
            // the interval properly intersect one another
            Relationship::Intersecting(interval) => {
                let mut maybe_start = None;
                for anti_diagonal_idx in interval.start..=interval.end {
                    let self_bound = self.get(anti_diagonal_idx);
                    let other_bound = other.get(anti_diagonal_idx);

                    if self_bound.intersects(other_bound) {
                        maybe_start = Some(anti_diagonal_idx);
                        break;
                    }
                }

                match maybe_start {
                    // if we've found the start of a proper intersection,
                    // we'll check to find the end of the intersection
                    Some(start) => {
                        let mut end = start;

                        for anti_diagonal_idx in (interval.start..=interval.end).rev() {
                            let self_bound = self.get(anti_diagonal_idx);
                            let other_bound = other.get(anti_diagonal_idx);

                            if self_bound.intersects(other_bound) {
                                end = anti_diagonal_idx;
                                break;
                            }
                        }
                        Relationship::Intersecting(Interval { start, end })
                    }
                    // if we have no start, then we also won't have an end
                    None => Relationship::Disjoint(Interval { start: 0, end: 0 }),
                }
            }
            disjoint => disjoint,
        }
    }

    pub fn ascii(&self, other: Option<&Self>) -> anyhow::Result<()> {
        let mut img = vec![vec![0; self.profile_length + 1]; self.target_length + 1];

        self.bounds().iter().for_each(|b| {
            b.cell_zip().for_each(|(t, p)| img[t][p] += 1);
        });

        if let Some(other) = other {
            if self.profile_length != other.profile_length
                || self.target_length != other.target_length
            {
                bail!("dimension mismatch")
            }

            other.bounds().iter().for_each(|b| {
                b.cell_zip().for_each(|(t, p)| img[t][p] += 1);
            });
        }

        print!("/// `     ");
        (0..=self.profile_length).for_each(|p| print!("{p}  "));
        println!("\x08`");
        println!("/// `    {} `", "---".repeat(self.profile_length + 1));
        for (row_idx, row) in img.iter().enumerate() {
            print!("/// ` {row_idx:2} |");
            for &v in row.iter() {
                if v == 0 {
                    print!("   ");
                } else if v == 1 {
                    print!("x  ");
                } else if v == 2 {
                    print!("*  ")
                } else {
                    bail!("cloud cell error")
                }
            }
            println!("`");
        }
        println!("/// `    {} `", "---".repeat(self.profile_length + 1));

        Ok(())
    }
}

#[cfg(feature = "debug")]
pub mod debug {
    use std::path::Path;

    use anyhow::Context;
    use image::{Rgb, RgbImage};

    use super::{AntiDiagonal, AntiDiagonalBounds};

    #[allow(dead_code)]
    pub fn cloud_image(
        bounds_a: &AntiDiagonalBounds,
        bounds_b: &AntiDiagonalBounds,
        path: impl AsRef<Path>,
    ) -> anyhow::Result<()> {
        let red = Rgb([255, 0, 0]);
        let blue = Rgb([0, 0, 255]);

        let width = bounds_a.profile_length;
        let height = bounds_a.target_length;

        let mut img = RgbImage::new(width as u32 + 1, height as u32 + 1);

        let mut draw = |bound: &AntiDiagonal, color: Rgb<u8>| {
            bound.cell_zip().for_each(|(y, x)| {
                img.put_pixel(x as u32, y as u32, color);
            });
        };

        bounds_a.bounds().iter().for_each(|b| draw(b, red));
        bounds_b.bounds().iter().for_each(|b| draw(b, blue));

        img.save(path).context("failed to write image")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use consts::*;

    use anyhow::anyhow;
    use pretty_assertions::assert_eq;

    mod consts {
        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |                                 `
        /// `  3 |         x  x  x  x              `
        /// `  4 |         x  x  x  x              `
        /// `  5 |         x  x  x  x              `
        /// `  6 |         x  x  x  x              `
        /// `  7 |                                 `
        /// `  8 |                                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_SQUARE_3_3_6_6: [[usize; 4]; 7] = [
            [3, 3, 3, 3],
            [4, 3, 3, 4],
            [5, 3, 3, 5],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [6, 5, 5, 6],
            [6, 6, 6, 6],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |            x                    `
        /// `  2 |      x  x  x  x  x              `
        /// `  3 |      x  x  x  x                 `
        /// `  4 |   x  x  x  x  x  x     x        `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |      x     x  x  x  x  x  x     `
        /// `  7 |               x  x  x  x        `
        /// `  8 |            x  x  x  x  x        `
        /// `  9 |                  x              `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_WINGS: [[usize; 4]; 13] = [
            [2, 2, 2, 2],
            [4, 1, 1, 4],
            [4, 2, 2, 4],
            [5, 2, 2, 5],
            [6, 2, 2, 6],
            [5, 4, 4, 5],
            [6, 4, 4, 6],
            [6, 5, 5, 6],
            [8, 4, 4, 8],
            [8, 5, 5, 8],
            [8, 6, 6, 8],
            [9, 6, 6, 9],
            [8, 8, 8, 8],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |      x  x  x  x                 `
        /// `  3 |      x  x  x  x                 `
        /// `  4 |      x  x  x  x  x              `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |            x  x  x  x  x        `
        /// `  7 |               x  x  x  x        `
        /// `  8 |               x  x  x  x        `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_WINGS_TRIMMED: [[usize; 4]; 13] = [
            [2, 2, 2, 2],
            [3, 2, 2, 3],
            [4, 2, 2, 4],
            [5, 2, 2, 5],
            [5, 3, 3, 5],
            [5, 4, 4, 5],
            [6, 4, 4, 6],
            [6, 5, 5, 6],
            [7, 5, 5, 7],
            [8, 5, 5, 8],
            [8, 6, 6, 8],
            [8, 7, 7, 8],
            [8, 8, 8, 8],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |                                 `
        /// `  3 |                  x              `
        /// `  4 |               x  x  x           `
        /// `  5 |            x  x  x              `
        /// `  6 |         x  x  x                 `
        /// `  7 |            x                    `
        /// `  8 |                                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_UNSQUARED: [[usize; 4]; 3] = [[6, 3, 3, 6], [6, 4, 4, 6], [7, 4, 4, 7]];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |                                 `
        /// `  3 |         x  x  x  x              `
        /// `  4 |         x  x  x  x  x           `
        /// `  5 |         x  x  x  x  x           `
        /// `  6 |         x  x  x  x  x           `
        /// `  7 |            x  x  x  x           `
        /// `  8 |                                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `

        pub const BOUNDS_SQUARED: [[usize; 4]; 9] = [
            [3, 3, 3, 3],
            [4, 3, 3, 4],
            [5, 3, 3, 5],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [7, 4, 4, 7],
            [7, 5, 5, 7],
            [7, 6, 6, 7],
            [7, 7, 7, 7],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x  x                    `
        /// `  2 |            x  x                 `
        /// `  3 |               x  x              `
        /// `  4 |                  x  x           `
        /// `  5 |                     x  x        `
        /// `  6 |                                 `
        /// `  7 |                                 `
        /// `  8 |                                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_A_1: [[usize; 4]; 12] = [
            [1, 1, 1, 1],
            [1, 2, 1, 2],
            [1, 3, 1, 3],
            [1, 4, 1, 4],
            [2, 4, 2, 4],
            [2, 5, 2, 5],
            [3, 5, 3, 5],
            [3, 6, 3, 6],
            [4, 6, 4, 6],
            [4, 7, 4, 7],
            [5, 7, 5, 7],
            [5, 8, 5, 8],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |                                 `
        /// `  3 |                                 `
        /// `  4 |                                 `
        /// `  5 |      x  x                       `
        /// `  6 |         x  x                    `
        /// `  7 |            x  x                 `
        /// `  8 |               x  x              `
        /// `  9 |                  x  x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_A_2: [[usize; 4]; 12] = [
            [5, 2, 5, 2],
            [5, 3, 5, 3],
            [6, 3, 6, 3],
            [6, 4, 6, 4],
            [7, 4, 7, 4],
            [7, 5, 7, 5],
            [8, 5, 8, 5],
            [8, 6, 8, 6],
            [9, 6, 9, 6],
            [9, 7, 9, 7],
            [9, 8, 9, 8],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x  x                    `
        /// `  2 |      x  x  x  x                 `
        /// `  3 |      x  x  x  x  x              `
        /// `  4 |      x  x  x  x  x  x           `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |         x  x  x  x  x  x        `
        /// `  7 |            x  x  x  x  x        `
        /// `  8 |               x  x  x  x        `
        /// `  9 |                  x  x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_A_MERGE: [[usize; 4]; 17] = [
            [1, 1, 1, 1],
            [1, 2, 1, 2],
            [2, 2, 1, 3],
            [3, 2, 1, 4],
            [4, 2, 2, 4],
            [5, 2, 2, 5],
            [5, 3, 3, 5],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [7, 4, 4, 7],
            [7, 5, 5, 7],
            [8, 5, 5, 8],
            [8, 6, 6, 8],
            [9, 6, 7, 8],
            [9, 7, 8, 8],
            [9, 8, 9, 8],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x                             `
        /// `  2 |   x                             `
        /// `  3 |   x                             `
        /// `  4 |   x                             `
        /// `  5 |   x  x                          `
        /// `  6 |      x  x                       `
        /// `  7 |         x  x                    `
        /// `  8 |            x  x                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_B_1: [[usize; 4]; 12] = [
            [1, 1, 1, 1],
            [2, 1, 2, 1],
            [3, 1, 3, 1],
            [4, 1, 4, 1],
            [5, 1, 5, 1],
            [5, 2, 5, 2],
            [6, 2, 6, 2],
            [6, 3, 6, 3],
            [7, 3, 7, 3],
            [7, 4, 7, 4],
            [8, 4, 8, 4],
            [8, 5, 8, 5],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |               x  x              `
        /// `  3 |                  x  x           `
        /// `  4 |                     x  x        `
        /// `  5 |                        x  x     `
        /// `  6 |                           x     `
        /// `  7 |                           x     `
        /// `  8 |                           x     `
        /// `  9 |                           x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_B_2: [[usize; 4]; 12] = [
            [2, 5, 2, 5],
            [2, 6, 2, 6],
            [3, 6, 3, 6],
            [3, 7, 3, 7],
            [4, 7, 4, 7],
            [4, 8, 4, 8],
            [5, 8, 5, 8],
            [5, 9, 5, 9],
            [6, 9, 6, 9],
            [7, 9, 7, 9],
            [8, 9, 8, 9],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x                             `
        /// `  2 |   x  x  x  x  x  x              `
        /// `  3 |   x  x  x  x  x  x  x           `
        /// `  4 |   x  x  x  x  x  x  x  x        `
        /// `  5 |   x  x  x  x  x  x  x  x  x     `
        /// `  6 |      x  x  x  x  x  x  x  x     `
        /// `  7 |         x  x  x  x  x  x  x     `
        /// `  8 |            x  x  x  x  x  x     `
        /// `  9 |                           x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_B_MERGE: [[usize; 4]; 17] = [
            [1, 1, 1, 1],
            [2, 1, 2, 1],
            [3, 1, 2, 2],
            [4, 1, 2, 3],
            [5, 1, 2, 4],
            [5, 2, 2, 5],
            [6, 2, 2, 6],
            [6, 3, 3, 6],
            [7, 3, 3, 7],
            [7, 4, 4, 7],
            [8, 4, 4, 8],
            [8, 5, 5, 8],
            [8, 6, 5, 9],
            [8, 7, 6, 9],
            [8, 8, 7, 9],
            [8, 9, 8, 9],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x                       `
        /// `  2 |   x  x  x  x  x                 `
        /// `  3 |   x  x  x  x  x  x              `
        /// `  4 |      x  x  x  x  x  x  x  x     `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |         x  x  x  x  x           `
        /// `  7 |            x  x  x              `
        /// `  8 |            x  x                 `
        /// `  9 |            x                    `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_C_1: [[usize; 4]; 12] = [
            [1, 1, 1, 1],
            [2, 1, 1, 2],
            [3, 1, 1, 3],
            [3, 2, 2, 3],
            [4, 2, 2, 4],
            [5, 2, 2, 5],
            [5, 3, 3, 5],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [7, 4, 4, 7],
            [8, 4, 4, 8],
            [9, 4, 4, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                  x              `
        /// `  2 |               x  x              `
        /// `  3 |            x  x  x              `
        /// `  4 |         x  x  x  x  x           `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |   x  x  x  x  x  x  x  x        `
        /// `  7 |            x  x  x  x  x  x     `
        /// `  8 |               x  x  x  x  x     `
        /// `  9 |                     x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_C_2: [[usize; 4]; 12] = [
            [6, 1, 1, 6],
            [6, 2, 2, 6],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [7, 4, 4, 7],
            [7, 5, 5, 7],
            [8, 5, 5, 8],
            [8, 6, 6, 8],
            [8, 7, 7, 8],
            [9, 7, 7, 9],
            [9, 8, 8, 9],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x  x  x  x              `
        /// `  2 |   x  x  x  x  x  x              `
        /// `  3 |   x  x  x  x  x  x              `
        /// `  4 |   x  x  x  x  x  x  x  x  x     `
        /// `  5 |   x  x  x  x  x  x  x  x  x     `
        /// `  6 |   x  x  x  x  x  x  x  x  x     `
        /// `  7 |            x  x  x  x  x  x     `
        /// `  8 |            x  x  x  x  x  x     `
        /// `  9 |            x  x  x  x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_C_MERGE: [[usize; 4]; 17] = [
            [1, 1, 1, 1],
            [2, 1, 1, 2],
            [3, 1, 1, 3],
            [4, 1, 1, 4],
            [5, 1, 1, 5],
            [6, 1, 1, 6],
            [6, 2, 2, 6],
            [6, 3, 3, 6],
            [6, 4, 4, 6],
            [7, 4, 4, 7],
            [8, 4, 4, 8],
            [9, 4, 4, 9],
            [9, 5, 5, 9],
            [9, 6, 6, 9],
            [9, 7, 7, 9],
            [9, 8, 8, 9],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x                       `
        /// `  2 |   x  x  x  x                    `
        /// `  3 |   x  x  x                       `
        /// `  4 |      x                          `
        /// `  5 |                                 `
        /// `  6 |                                 `
        /// `  7 |                                 `
        /// `  8 |                                 `
        /// `  9 |                                 `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_D_1: [[usize; 4]; 5] = [
            [1, 1, 1, 1],
            [2, 1, 1, 2],
            [3, 1, 1, 3],
            [3, 2, 2, 3],
            [4, 2, 2, 4],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |                                 `
        /// `  2 |                                 `
        /// `  3 |                                 `
        /// `  4 |                                 `
        /// `  5 |                                 `
        /// `  6 |                        x        `
        /// `  7 |                     x  x  x     `
        /// `  8 |                  x  x  x  x     `
        /// `  9 |                     x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_D_2: [[usize; 4]; 5] = [
            [8, 6, 6, 8],
            [8, 7, 7, 8],
            [9, 7, 7, 9],
            [9, 8, 8, 9],
            [9, 9, 9, 9],
        ];

        /// `     0  1  2  3  4  5  6  7  8  9  10 `
        /// `    --------------------------------- `
        /// `  0 |                                 `
        /// `  1 |   x  x  x                       `
        /// `  2 |   x  x  x  x  x  x  x  x        `
        /// `  3 |   x  x  x  x  x  x  x  x        `
        /// `  4 |      x  x  x  x  x  x  x        `
        /// `  5 |      x  x  x  x  x  x  x        `
        /// `  6 |      x  x  x  x  x  x  x        `
        /// `  7 |      x  x  x  x  x  x  x  x     `
        /// `  8 |      x  x  x  x  x  x  x  x     `
        /// `  9 |                     x  x  x     `
        /// ` 10 |                                 `
        /// `    --------------------------------- `
        pub const BOUNDS_D_MERGE: [[usize; 4]; 17] = [
            [1, 1, 1, 1],
            [2, 1, 1, 2],
            [3, 1, 1, 3],
            [3, 2, 2, 3],
            [4, 2, 2, 4],
            [5, 2, 2, 5],
            [6, 2, 2, 6],
            [7, 2, 2, 7],
            [8, 2, 2, 8],
            [8, 3, 3, 8],
            [8, 4, 4, 8],
            [8, 5, 5, 8],
            [8, 6, 6, 8],
            [8, 7, 7, 8],
            [9, 7, 7, 9],
            [9, 8, 8, 9],
            [9, 9, 9, 9],
        ];
    }

    /// This is a helper function for hand-filling AD bounds.
    fn fill_bounds(bounds: &mut AntiDiagonalBounds, slice: &[[usize; 4]]) -> anyhow::Result<()> {
        bounds.reset();

        let first = &slice.first().ok_or(anyhow!("empty vec"))?;
        let last = &slice.last().ok_or(anyhow!("empty vec"))?;

        let first_idx = first[0] + first[1];
        let last_idx = last[0] + last[1];

        // this will store the number of times each anti-diagonal is set
        let mut counts = vec![0; last_idx + 1];

        for e in slice.iter() {
            let a = e[0] + e[1];
            let b = e[2] + e[3];

            if a != b {
                bail!("index mismatch: {a} != {b}")
            } else if a < first_idx || a > last_idx || b < first_idx || b > last_idx {
                bail!("index out of range: {e:?} {first_idx} {last_idx}")
            }

            counts[a] += 1;
        }

        for (idx, &count) in counts.iter().enumerate() {
            let in_interval = idx >= first_idx && idx <= last_idx;

            if (!in_interval && count != 0) || (in_interval && count != 1) {
                bail!("index error: AD:{idx} ({count})")
            }
        }

        bounds.min_anti_diagonal_idx = first_idx;
        bounds.max_anti_diagonal_idx = last_idx;

        slice.iter().for_each(|e| {
            let idx = e[0] + e[1];
            bounds.set(idx, e[0], e[1], e[2], e[3]);
        });

        bounds.bounds().iter().try_for_each(|b| {
            if !b.valid() {
                Err(anyhow!("invalid bound: {b:?}"))
            } else {
                Ok(())
            }
        })
    }

    #[test]
    pub fn test_advance_forward() {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.set(2, 1, 1, 1, 1);
        bounds.min_anti_diagonal_idx = 2;
        bounds.max_anti_diagonal_idx = 2;

        (0..18).for_each(|_| bounds.advance_forward());

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        target_bounds.fill_rectangle(1, 1, 10, 10);

        assert_eq!(bounds, target_bounds);
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

        assert_eq!(bounds, target_bounds);
    }

    #[test]
    pub fn test_fill_square() -> anyhow::Result<()> {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        bounds.fill_rectangle(3, 3, 6, 6);

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);

        fill_bounds(&mut target_bounds, &BOUNDS_SQUARE_3_3_6_6)?;
        assert_eq!(bounds, target_bounds);

        Ok(())
    }

    #[test]
    pub fn test_trim_wings() -> anyhow::Result<()> {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        fill_bounds(&mut bounds, &BOUNDS_WINGS)?;

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        fill_bounds(&mut target_bounds, &BOUNDS_WINGS_TRIMMED)?;

        bounds.trim_wings()?;
        assert_eq!(bounds, target_bounds);

        // re-run it to make sure it does
        // nothing to already-trimmed bounds
        bounds.trim_wings()?;
        assert_eq!(bounds, target_bounds);

        Ok(())
    }

    #[test]
    fn test_square_corners() -> anyhow::Result<()> {
        let mut bounds = AntiDiagonalBounds::new(10, 10);
        fill_bounds(&mut bounds, &BOUNDS_UNSQUARED)?;

        let mut target_bounds = AntiDiagonalBounds::new(10, 10);
        fill_bounds(&mut target_bounds, &BOUNDS_SQUARED)?;

        target_bounds.min_anti_diagonal_idx = 6;
        target_bounds.max_anti_diagonal_idx = 14;

        bounds.square_corners();
        assert_eq!(bounds, target_bounds);

        // re-run it to make sure it does
        // nothing to already-squared bounds
        bounds.square_corners();
        assert_eq!(bounds, target_bounds);

        Ok(())
    }

    #[test]
    fn test_merge() -> anyhow::Result<()> {
        let mut b1 = AntiDiagonalBounds::new(10, 10);
        let mut b2 = AntiDiagonalBounds::new(10, 10);
        let mut target_bounds = AntiDiagonalBounds::new(10, 10);

        //--------------------------------
        fill_bounds(&mut b1, &BOUNDS_A_1)?;
        fill_bounds(&mut b2, &BOUNDS_A_2)?;
        fill_bounds(&mut target_bounds, &BOUNDS_A_MERGE)?;

        b1.merge(&b2);
        assert_eq!(b1, target_bounds);

        // reset and invert merge order
        fill_bounds(&mut b1, &BOUNDS_A_1)?;

        b2.merge(&b1);
        assert_eq!(b2, target_bounds);

        //--------------------------------
        fill_bounds(&mut b1, &BOUNDS_B_1)?;
        fill_bounds(&mut b2, &BOUNDS_B_2)?;
        fill_bounds(&mut target_bounds, &BOUNDS_B_MERGE)?;

        b1.merge(&b2);
        assert_eq!(b1, target_bounds);

        // reset and invert merge order
        fill_bounds(&mut b1, &BOUNDS_B_1)?;

        b2.merge(&b1);
        assert_eq!(b2, target_bounds);

        //--------------------------------
        fill_bounds(&mut b1, &BOUNDS_C_1)?;
        fill_bounds(&mut b2, &BOUNDS_C_2)?;
        fill_bounds(&mut target_bounds, &BOUNDS_C_MERGE)?;

        b1.merge(&b2);
        assert_eq!(b1, target_bounds);

        // reset and invert merge order
        fill_bounds(&mut b1, &BOUNDS_C_1)?;

        b2.merge(&b1);
        assert_eq!(b2, target_bounds);

        //--------------------------------
        fill_bounds(&mut b1, &BOUNDS_D_1)?;
        fill_bounds(&mut b2, &BOUNDS_D_2)?;
        fill_bounds(&mut target_bounds, &BOUNDS_D_MERGE)?;

        b1.ascii(None);
        b2.ascii(None);

        b1.merge(&b2);
        b1.ascii(None);

        assert_eq!(b1, target_bounds);

        // reset and invert merge order
        // fill_bounds(&mut b1, &BOUNDS_D_1)?;

        b2.merge(&b1);
        assert_eq!(b2, target_bounds);
        Ok(())
    }

    #[test]
    fn test_anti_diagonal_grow_up() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_up(6);
        assert_eq!(a.right_target_idx, 5);
        assert_eq!(a.right_profile_idx, 5);

        // make sure we grow properly
        a.grow_up(1);
        assert_eq!(a.right_target_idx, 1);
        assert_eq!(a.right_profile_idx, 9);

        // make sure we can't grow past 1
        a.grow_up(0);
        assert_eq!(a.right_target_idx, 1);
        assert_eq!(a.right_profile_idx, 9);
    }

    #[test]
    fn test_anti_diagonal_grow_down() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_down(4);
        assert_eq!(a.left_target_idx, 5);
        assert_eq!(a.left_profile_idx, 5);

        // make sure we grow properly
        a.grow_down(9);
        assert_eq!(a.left_target_idx, 9);
        assert_eq!(a.left_profile_idx, 1);
    }

    #[test]
    fn test_anti_diagonal_grow_left() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_left(6);
        assert_eq!(a.left_target_idx, 5);
        assert_eq!(a.left_profile_idx, 5);

        // make sure we grow properly
        a.grow_left(1);
        assert_eq!(a.left_target_idx, 9);
        assert_eq!(a.left_profile_idx, 1);

        // make sure we can't grow past 1
        a.grow_left(0);
        assert_eq!(a.left_target_idx, 9);
        assert_eq!(a.left_profile_idx, 1);
    }

    #[test]
    fn test_anti_diagonal_grow_right() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_right(4);
        assert_eq!(a.right_target_idx, 5);
        assert_eq!(a.right_profile_idx, 5);

        // make sure we grow properly
        a.grow_right(9);
        assert_eq!(a.right_target_idx, 1);
        assert_eq!(a.right_profile_idx, 9);
    }

    #[test]
    fn test_anti_diagonal_intersects() {
        // different anti-diagonals
        let a = AntiDiagonal::new(9, 1, 1, 9);
        let b = AntiDiagonal::new(9, 2, 2, 9);
        assert!(!a.intersects(&b));

        // same antidiagonal ovlerap
        let a = AntiDiagonal::new(9, 1, 1, 9);
        let b = AntiDiagonal::new(9, 1, 1, 9);
        assert!(a.intersects(&b));

        // single cell overlap
        let a = AntiDiagonal::new(9, 1, 5, 5);
        let b = AntiDiagonal::new(5, 5, 1, 9);
        assert!(a.intersects(&b));

        // single cell difference
        let a = AntiDiagonal::new(9, 1, 6, 4);
        let b = AntiDiagonal::new(5, 5, 1, 9);
        assert!(!a.intersects(&b));

        // several cell overlap
        let a = AntiDiagonal::new(9, 1, 3, 7);
        let b = AntiDiagonal::new(7, 3, 1, 9);
        assert!(a.intersects(&b));

        // several cell difference
        let a = AntiDiagonal::new(9, 1, 7, 3);
        let b = AntiDiagonal::new(3, 7, 1, 9);
        assert!(!a.intersects(&b));
    }

    #[test]
    fn test_cloud_relationship() {
        // no overlap, not touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 4, 4);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(6, 6, 10, 10);

        // bounds should always intersect with themselves
        assert_eq!(
            a.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 2, end: 8 })
        );
        assert_eq!(
            b.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 12, end: 20 })
        );

        assert_eq!(
            a.cloud_relationship(&b),
            Relationship::Disjoint(Interval { start: 9, end: 11 })
        );
        assert_eq!(
            b.cloud_relationship(&a),
            Relationship::Disjoint(Interval { start: 9, end: 11 })
        );

        // --------------------------------------------------------
        // no overlap, touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 5, 6);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(6, 5, 10, 10);

        // bounds should always intersect with themselves
        assert_eq!(
            a.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 2, end: 11 })
        );
        assert_eq!(
            b.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 11, end: 20 })
        );

        assert_eq!(
            a.cloud_relationship(&b),
            Relationship::Disjoint(Interval { start: 0, end: 0 }),
        );
        assert_eq!(
            b.cloud_relationship(&a),
            Relationship::Disjoint(Interval { start: 0, end: 0 }),
        );

        // --------------------------------------------------------
        // 1 cell overlap
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 5, 5);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(5, 5, 10, 10);

        // bounds should always intersect with themselves
        assert_eq!(
            a.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 2, end: 10 }),
        );
        assert_eq!(
            b.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 10, end: 20 }),
        );

        assert_eq!(
            a.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 10, end: 10 }),
        );
        assert_eq!(
            b.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 10, end: 10 }),
        );

        // --------------------------------------------------------
        // many cell overlap
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 8, 8);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(2, 2, 10, 10);

        // bounds should always intersect with themselves
        assert_eq!(
            a.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 2, end: 16 }),
        );
        assert_eq!(
            b.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 4, end: 20 }),
        );

        assert_eq!(
            a.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 4, end: 16 }),
        );
        assert_eq!(
            b.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 4, end: 16 }),
        );

        // --------------------------------------------------------
        // anti-diagonal overlap, not touching
        let mut a = AntiDiagonalBounds::new(10, 10);
        a.fill_rectangle(1, 1, 8, 4);

        let mut b = AntiDiagonalBounds::new(10, 10);
        b.fill_rectangle(1, 5, 8, 9);

        // bounds should always intersect with themselves
        assert_eq!(
            a.cloud_relationship(&a),
            Relationship::Intersecting(Interval { start: 2, end: 12 }),
        );
        assert_eq!(
            b.cloud_relationship(&b),
            Relationship::Intersecting(Interval { start: 6, end: 17 }),
        );

        assert_eq!(
            a.cloud_relationship(&b),
            Relationship::Disjoint(Interval { start: 0, end: 0 }),
        );
        assert_eq!(
            b.cloud_relationship(&a),
            Relationship::Disjoint(Interval { start: 0, end: 0 }),
        );
    }
}
