use anyhow::{anyhow, bail};

use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::io::{stdout, Write};
use std::iter::{Rev, Zip};
use std::ops::RangeInclusive;

use super::{BackgroundCell, BackgroundState, CoreCell, CoreState};

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub fn length(&self) -> usize {
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

#[derive(Default, Clone)]
pub struct Cell {
    pub profile_idx: usize,
    pub seq_idx: usize,
}

impl Cell {
    pub fn idx(&self) -> usize {
        self.seq_idx + self.profile_idx
    }

    pub fn m_cell(&self) -> CoreCell {
        (CoreState::M(self.profile_idx), self.seq_idx)
    }

    pub fn i_cell(&self) -> CoreCell {
        (CoreState::I(self.profile_idx), self.seq_idx)
    }

    pub fn d_cell(&self) -> CoreCell {
        (CoreState::D(self.profile_idx), self.seq_idx)
    }

    pub fn b_cell(&self) -> BackgroundCell {
        (BackgroundState::B, self.seq_idx)
    }

    pub fn e_cell(&self) -> BackgroundCell {
        (BackgroundState::E, self.seq_idx)
    }
}

impl Debug for Cell {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "(p{}, s{})", self.profile_idx, self.seq_idx)
    }
}

#[derive(Default, Clone)]
pub struct Bound(Cell, Cell);
pub struct ArrayProfileMajorBound([usize; 4]);
pub struct ArraySeqMajorBound([usize; 4]);
pub struct TupleProfileMajorBound((usize, usize), (usize, usize));
pub struct TupleSeqMajorBound((usize, usize), (usize, usize));

pub struct BoundIter {
    profile_start: usize,
    profile_end: usize,
    seq_start: usize,
    seq_end: usize,
}

impl Iterator for BoundIter {
    type Item = Cell;

    fn next(&mut self) -> Option<Self::Item> {
        if self.profile_start < self.profile_end {
            let cell = Cell {
                profile_idx: self.profile_start,
                seq_idx: self.seq_start,
            };
            self.profile_start -= 1;
            self.seq_start += 1;
            Some(cell)
        } else {
            None
        }
    }
}

impl DoubleEndedIterator for BoundIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.profile_start < self.profile_end {
            let cell = Cell {
                profile_idx: self.profile_end,
                seq_idx: self.seq_end,
            };
            self.profile_start += 1;
            self.seq_start -= 1;
            Some(cell)
        } else {
            None
        }
    }
}

impl Bound {
    pub fn idx(&self) -> usize {
        debug_assert!(self.0.idx() == self.1.idx());
        self.0.idx()
    }

    pub fn iter(&self) -> BoundIter {
        BoundIter {
            profile_start: self.1.profile_idx,
            profile_end: self.0.profile_idx,
            seq_start: self.0.seq_idx,
            seq_end: self.1.seq_idx,
        }
    }

    pub fn from_anti_diagonal(ad: &AntiDiagonal) -> Self {
        Self(
            Cell {
                profile_idx: ad.right_profile_idx,
                seq_idx: ad.right_target_idx,
            },
            Cell {
                profile_idx: ad.left_profile_idx,
                seq_idx: ad.left_target_idx,
            },
        )
    }
}

impl Debug for Bound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(p{}, s{}), (p{}, s{}) | AD: {}",
            self.0.profile_idx,
            self.0.seq_idx,
            self.1.profile_idx,
            self.1.seq_idx,
            self.idx()
        )
    }
}

pub trait BoundInterpretable {
    fn interpret(&self) -> Bound;
}

impl BoundInterpretable for Bound {
    fn interpret(&self) -> Bound {
        self.clone()
    }
}

impl BoundInterpretable for ArrayProfileMajorBound {
    fn interpret(&self) -> Bound {
        Bound(
            Cell {
                profile_idx: self.0[0],
                seq_idx: self.0[1],
            },
            Cell {
                profile_idx: self.0[2],
                seq_idx: self.0[3],
            },
        )
    }
}

impl BoundInterpretable for ArraySeqMajorBound {
    fn interpret(&self) -> Bound {
        Bound(
            Cell {
                profile_idx: self.0[3],
                seq_idx: self.0[2],
            },
            Cell {
                profile_idx: self.0[1],
                seq_idx: self.0[0],
            },
        )
    }
}

impl BoundInterpretable for TupleProfileMajorBound {
    fn interpret(&self) -> Bound {
        Bound(
            Cell {
                profile_idx: self.0 .0,
                seq_idx: self.0 .1,
            },
            Cell {
                profile_idx: self.1 .0,
                seq_idx: self.1 .1,
            },
        )
    }
}

impl BoundInterpretable for TupleSeqMajorBound {
    fn interpret(&self) -> Bound {
        Bound(
            Cell {
                profile_idx: self.1 .1,
                seq_idx: self.1 .0,
            },
            Cell {
                profile_idx: self.0 .1,
                seq_idx: self.0 .0,
            },
        )
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

    // TODO: better API for this
    pub fn left_target_major(&self) -> Cell {
        Cell {
            seq_idx: self.left_target_idx,
            profile_idx: self.left_profile_idx,
        }
    }

    pub fn right_target_major(&self) -> Cell {
        Cell {
            seq_idx: self.right_target_idx,
            profile_idx: self.right_profile_idx,
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

    pub fn len(&self) -> usize {
        self.right_profile_idx - self.left_profile_idx + 1
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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

    pub fn grow_target(&mut self, target_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.right_target_idx = self.right_target_idx.min(target_idx).max(1);
        self.right_profile_idx = idx - self.right_target_idx;

        debug_assert!(idx == self.idx())
    }

    pub fn grow_down(&mut self, target_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.left_target_idx = self.left_target_idx.max(target_idx);
        self.left_profile_idx = idx - self.left_target_idx;

        debug_assert!(idx == self.idx())
    }

    pub fn grow_left(&mut self, profile_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.left_profile_idx = self.left_profile_idx.min(profile_idx).max(1);
        self.left_target_idx = idx - self.left_profile_idx;

        debug_assert!(idx == self.idx())
    }

    pub fn grow_right(&mut self, profile_idx: usize) {
        debug_assert!(self.valid());

        let idx = self.idx();
        self.right_profile_idx = self.right_profile_idx.max(profile_idx);
        self.right_target_idx = idx - self.right_profile_idx;

        debug_assert!(idx == self.idx())
    }
}

pub enum CloudValidity {
    Valid,
    BoundInvalid {
        bound: Bound,
    },
    BoundOutOfRange {
        ad_start: usize,
        ad_end: usize,
        bound: Bound,
    },
}

#[derive(Clone, PartialEq, Debug)]
pub struct Cloud {
    pub bounds: Vec<AntiDiagonal>,
    pub target_length: usize,
    pub profile_length: usize,
    pub size: usize,
    pub ad_start: usize,
    pub ad_end: usize,
}

impl Default for Cloud {
    fn default() -> Self {
        Self {
            bounds: vec![AntiDiagonal::default()],
            target_length: 0,
            profile_length: 0,
            size: 1,
            ad_start: 0,
            ad_end: 0,
        }
    }
}

impl Cloud {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        let size = target_length + profile_length + 1;
        Self {
            bounds: vec![AntiDiagonal::default(); size],
            target_length,
            profile_length,
            size,
            ad_start: 0,
            ad_end: 0,
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

        self.ad_start = new_size;
        self.ad_end = 0;

        self.target_length = target_length;
        self.profile_length = profile_length;
    }

    /// Fill a Cloud with a list of hand-written bounds.
    ///
    /// For now, this is not part of the public API and is mainly
    /// a helper function for writing tests involving Clouds.
    #[allow(dead_code)]
    pub(crate) fn fill(&mut self, bounds: &[impl BoundInterpretable]) -> anyhow::Result<()> {
        let bounds: Vec<Bound> = bounds.iter().map(|b| b.interpret()).collect();
        self.reset();

        let first = &bounds.first().ok_or(anyhow!("empty vec"))?;
        let last = &bounds.last().ok_or(anyhow!("empty vec"))?;

        self.ad_start = first.idx();
        self.ad_end = last.idx();

        bounds.iter().for_each(|bound| {
            self.set(
                bound.idx(),
                bound.1.seq_idx,
                bound.1.profile_idx,
                bound.0.seq_idx,
                bound.0.profile_idx,
            );
        });

        match self.validity() {
            CloudValidity::Valid => Ok(()),
            CloudValidity::BoundInvalid { bound } => Err(anyhow!("invalid bound: {bound:?}")),
            CloudValidity::BoundOutOfRange {
                ad_start,
                ad_end,
                bound,
            } => Err(anyhow!(
                "bound: {bound:?}\nout of range: [{ad_start}, {ad_end}]",
            )),
        }
    }

    pub fn target_len(&self) -> usize {
        self.target_length
    }

    pub fn profile_len(&self) -> usize {
        self.profile_length
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
        // debug_assert!(
        //     left_target_idx + left_profile_idx,
        //     right_target_idx + right_profile_idx
        // );

        // make sure we are setting the anti-diagonal that we think we are
        //debug_assert!(anti_diagonal_idx, left_target_idx + left_profile_idx);

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
        self.get(self.ad_start)
    }

    pub fn last(&self) -> &AntiDiagonal {
        self.get(self.ad_end)
    }

    pub fn valid(&self) -> bool {
        self.bounds().iter().all(|bound| bound.valid())
    }

    pub fn validity(&self) -> CloudValidity {
        let mut out_of_range = self.bounds[0..self.ad_start]
            .iter()
            .chain(self.bounds[self.ad_end + 1..].iter());

        let mut in_range = self.bounds().iter();

        if let Some(ad) = out_of_range.find(|b| !b.is_default()) {
            CloudValidity::BoundOutOfRange {
                ad_start: self.ad_start,
                ad_end: self.ad_end,
                bound: Bound::from_anti_diagonal(ad),
            }
        } else if let Some(ad) = in_range.find(|b| !b.valid()) {
            CloudValidity::BoundInvalid {
                bound: Bound::from_anti_diagonal(ad),
            }
        } else {
            CloudValidity::Valid
        }
    }

    /// Get the total number of cells that exist within the cloud boundaries.
    pub fn cloud_size(&self) -> usize {
        let mut cloud_size = 0usize;
        for bound in self.bounds[self.ad_start..=self.ad_end].iter() {
            cloud_size += bound.left_target_idx - bound.right_target_idx;
        }
        cloud_size
    }

    /// Get the number of anti-diagonals defined in the cloud.
    pub fn num_anti_diagonals(&self) -> usize {
        self.ad_end - self.ad_start + 1
    }

    pub fn bounds(&self) -> &[AntiDiagonal] {
        &self.bounds[self.ad_start..=self.ad_end]
    }

    pub fn bounds_mut(&mut self) -> &mut [AntiDiagonal] {
        &mut self.bounds[self.ad_start..=self.ad_end]
    }

    /// This removes all of the protruding regions in the cloud that are
    /// unreachable from a traceback that traverses the entire cloud.
    pub fn trim_wings(&mut self) -> anyhow::Result<()> {
        for anti_diagonal_idx in self.ad_start + 1..=self.ad_end {
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

        for anti_diagonal_idx in (self.ad_start..self.ad_end).rev() {
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

        self.ad_start = first_anti_diagonal_idx - left_distance;
        self.ad_end = last_anti_diagonal_idx + right_distance;
    }

    pub fn advance_forward(&mut self) {
        let last_bound = self.last();

        let next_idx = self.ad_end + 1;

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

        self.ad_end = next_idx;
    }

    pub fn advance_reverse(&mut self) {
        let first_bound = self.first();

        let next_idx = self.ad_start - 1;

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

        self.ad_start = next_idx;
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

        self.ad_start = anti_diagonal_start;
        self.ad_end = anti_diagonal_end;
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
        self.ad_start.cmp(&other.ad_start)
    }

    pub fn right_offset_to(&mut self, other: &Self) -> Ordering {
        self.ad_end.cmp(&other.ad_end)
    }

    pub fn merge(&mut self, other: &Self) {
        let relationship = self.anti_diagonal_relationship(other);
        let interval = relationship.interval();

        // if 'self' is right-offset from 'other', we need
        // to pull in the remaining 'other' bounds
        if let Ordering::Greater = self.left_offset_to(other) {
            let self_slice = &mut self.bounds[other.ad_start..interval.start];
            let other_slice = &other.bounds[other.ad_start..interval.start];

            self_slice
                .iter_mut()
                .zip(other_slice)
                .for_each(|(self_bound, other_bound)| {
                    self_bound.replace_with(other_bound);
                });

            self.ad_start = other.ad_start;
        }

        // if 'self' is left-offset from 'other', we need
        // to pull in the remaining 'other' bounds
        if let Ordering::Less = self.right_offset_to(other) {
            let self_slice = &mut self.bounds[(interval.end + 1)..=other.ad_end];
            let other_slice = &other.bounds[(interval.end + 1)..=other.ad_end];

            self_slice
                .iter_mut()
                .zip(other_slice)
                .for_each(|(self_bound, other_bound)| {
                    self_bound.replace_with(other_bound);
                });

            self.ad_end = other.ad_end;
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

                let self_slice = &mut self.bounds[self.ad_start..interval.start];

                self_slice.iter_mut().for_each(|bound| {
                    bound.grow_target(min_target_idx);
                    bound.grow_left(min_profile_idx);

                    debug_assert!(bound.right_target_idx >= 1);
                    debug_assert!(bound.left_profile_idx >= 1);
                });

                let max_target_idx = self.get(interval.end).left_target_idx;
                let max_profile_idx = self.get(interval.end).right_profile_idx;

                let self_slice = &mut self.bounds[(interval.end + 1)..=self.ad_end];

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
                    .take(interval.length())
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
        let overlap_start = self.ad_start.max(other.ad_start);
        let overlap_end = self.ad_end.min(other.ad_end);

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

    pub fn vec_image(&self, other: Option<&Self>) -> anyhow::Result<Vec<Vec<u8>>> {
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
                b.cell_zip().for_each(|(t, p)| img[t][p] += 2);
            });
        };
        Ok(img)
    }

    #[allow(dead_code)]
    pub(crate) fn image(&self) -> CloudImage {
        let img = self.vec_image(None).unwrap();
        CloudImage {
            img_bytes: (0..img.iter().map(Vec::len).max().unwrap_or(0))
                .map(|i| img.iter().filter_map(|row| row.get(i).copied()).collect())
                .collect(),
            print_orientation: ImageOrientation::ProfileMajor,
            print_mode: PrintMode::Default,
        }
    }
}

pub enum ImageOrientation {
    ProfileMajor,
    SeqMajor,
}

pub enum PrintMode {
    Default,
    DocString,
}

pub(crate) struct CloudImage {
    img_bytes: Vec<Vec<u8>>,
    print_orientation: ImageOrientation,
    print_mode: PrintMode,
}

#[allow(dead_code)]
impl CloudImage {
    const BACKSPACE: &str = "\x08";

    pub fn profile_len(&self) -> usize {
        self.img_bytes.len()
    }

    pub fn seq_len(&self) -> usize {
        self.img_bytes[0].len()
    }

    pub fn profile_major(mut self) -> Self {
        self.print_orientation = ImageOrientation::ProfileMajor;
        self
    }

    pub fn seq_major(mut self) -> Self {
        self.print_orientation = ImageOrientation::SeqMajor;
        self
    }

    pub fn doc_string(mut self) -> Self {
        self.print_mode = PrintMode::DocString;
        self
    }

    pub fn print(&self) {
        self.write(&mut stdout()).unwrap()
    }

    fn write(&self, buf: &mut impl Write) -> anyhow::Result<()> {
        let (prefix, suffix) = match self.print_mode {
            PrintMode::Default => ("", ""),
            PrintMode::DocString => ("/// `", "`"),
        };

        let (rows, n_cols) = match self.print_orientation {
            ImageOrientation::ProfileMajor => (self.img_bytes.clone(), self.profile_len()),
            ImageOrientation::SeqMajor => (
                (0..self.img_bytes.iter().map(Vec::len).max().unwrap_or(0))
                    .map(|i| {
                        self.img_bytes
                            .iter()
                            .filter_map(|row| row.get(i).copied())
                            .collect()
                    })
                    .collect(),
                self.seq_len(),
            ),
        };

        // write the column indices
        write!(buf, "{}     ", prefix)?;
        (0..n_cols).try_for_each(|col_idx| write!(buf, "{col_idx}  "))?;
        println!("{}{}", Self::BACKSPACE, suffix);
        println!("{}    {} {}", prefix, "---".repeat(n_cols), suffix);

        // write the data
        for (row_idx, row) in rows.iter().enumerate() {
            write!(buf, "{} {row_idx:2} |", prefix)?;
            for &v in row.iter() {
                match v {
                    0 => write!(buf, "-  "),
                    1 => write!(buf, "x  "),
                    2 => write!(buf, "*  "),
                    3 => write!(buf, "o  "),
                    _ => bail!("cloud cell error"),
                }?
            }
            println!("{suffix}");
        }

        println!("{}    {} {}", prefix, "---".repeat(n_cols), suffix);
        Ok(())
    }
}

#[cfg(feature = "debug")]
pub mod debug {
    use std::path::Path;

    use anyhow::{bail, Context};
    use image::RgbImage;

    use super::Cloud;
    use consts::*;

    mod consts {
        use image::Rgb;

        pub const RED: Rgb<u8> = Rgb([255, 0, 0]);
        pub const GREEN: Rgb<u8> = Rgb([0, 255, 0]);
        pub const BLUE: Rgb<u8> = Rgb([0, 0, 255]);
    }

    impl Cloud {
        #[allow(dead_code)]
        pub fn write_image(
            &self,
            other: Option<&Self>,
            path: impl AsRef<Path>,
        ) -> anyhow::Result<()> {
            let vec_img = self.vec_image(other)?;

            let mut img = RgbImage::new(vec_img[0].len() as u32, vec_img.len() as u32);

            for (y, row) in vec_img.iter().enumerate() {
                for (x, &val) in row.iter().enumerate() {
                    if val == 0 {
                        // nothing
                    } else if val == 1 {
                        img.put_pixel(x as u32, y as u32, BLUE)
                    } else if val == 2 {
                        img.put_pixel(x as u32, y as u32, RED)
                    } else if val == 3 {
                        img.put_pixel(x as u32, y as u32, GREEN)
                    } else {
                        bail!("cloud cell error")
                    }
                }
            }

            img.save(path).context("failed to write image")
        }
    }
}

#[cfg(test)]
pub(crate) mod test_consts {
    use super::ArraySeqMajorBound;

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
    pub const BOUNDS_SQUARE_3_3_6_6: [ArraySeqMajorBound; 7] = [
        ArraySeqMajorBound([3, 3, 3, 3]),
        ArraySeqMajorBound([4, 3, 3, 4]),
        ArraySeqMajorBound([5, 3, 3, 5]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([6, 5, 5, 6]),
        ArraySeqMajorBound([6, 6, 6, 6]),
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
    pub const BOUNDS_WINGS: [ArraySeqMajorBound; 13] = [
        ArraySeqMajorBound([2, 2, 2, 2]),
        ArraySeqMajorBound([4, 1, 1, 4]),
        ArraySeqMajorBound([4, 2, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([6, 2, 2, 6]),
        ArraySeqMajorBound([5, 4, 4, 5]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([6, 5, 5, 6]),
        ArraySeqMajorBound([8, 4, 4, 8]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([9, 6, 6, 9]),
        ArraySeqMajorBound([8, 8, 8, 8]),
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
    pub const BOUNDS_WINGS_TRIMMED: [ArraySeqMajorBound; 13] = [
        ArraySeqMajorBound([2, 2, 2, 2]),
        ArraySeqMajorBound([3, 2, 2, 3]),
        ArraySeqMajorBound([4, 2, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([5, 3, 3, 5]),
        ArraySeqMajorBound([5, 4, 4, 5]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([6, 5, 5, 6]),
        ArraySeqMajorBound([7, 5, 5, 7]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([8, 7, 7, 8]),
        ArraySeqMajorBound([8, 8, 8, 8]),
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
    pub const BOUNDS_UNSQUARED: [ArraySeqMajorBound; 3] = [
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
    ];

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
    pub const BOUNDS_SQUARED: [ArraySeqMajorBound; 9] = [
        ArraySeqMajorBound([3, 3, 3, 3]),
        ArraySeqMajorBound([4, 3, 3, 4]),
        ArraySeqMajorBound([5, 3, 3, 5]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([7, 5, 5, 7]),
        ArraySeqMajorBound([7, 6, 6, 7]),
        ArraySeqMajorBound([7, 7, 7, 7]),
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
    pub const BOUNDS_A_1: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([1, 2, 1, 2]),
        ArraySeqMajorBound([1, 3, 1, 3]),
        ArraySeqMajorBound([1, 4, 1, 4]),
        ArraySeqMajorBound([2, 4, 2, 4]),
        ArraySeqMajorBound([2, 5, 2, 5]),
        ArraySeqMajorBound([3, 5, 3, 5]),
        ArraySeqMajorBound([3, 6, 3, 6]),
        ArraySeqMajorBound([4, 6, 4, 6]),
        ArraySeqMajorBound([4, 7, 4, 7]),
        ArraySeqMajorBound([5, 7, 5, 7]),
        ArraySeqMajorBound([5, 8, 5, 8]),
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
    pub const BOUNDS_A_2: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([5, 2, 5, 2]),
        ArraySeqMajorBound([5, 3, 5, 3]),
        ArraySeqMajorBound([6, 3, 6, 3]),
        ArraySeqMajorBound([6, 4, 6, 4]),
        ArraySeqMajorBound([7, 4, 7, 4]),
        ArraySeqMajorBound([7, 5, 7, 5]),
        ArraySeqMajorBound([8, 5, 8, 5]),
        ArraySeqMajorBound([8, 6, 8, 6]),
        ArraySeqMajorBound([9, 6, 9, 6]),
        ArraySeqMajorBound([9, 7, 9, 7]),
        ArraySeqMajorBound([9, 8, 9, 8]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_A_MERGE: [ArraySeqMajorBound; 17] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([1, 2, 1, 2]),
        ArraySeqMajorBound([2, 2, 1, 3]),
        ArraySeqMajorBound([3, 2, 1, 4]),
        ArraySeqMajorBound([4, 2, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([5, 3, 3, 5]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([7, 5, 5, 7]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([9, 6, 7, 8]),
        ArraySeqMajorBound([9, 7, 8, 8]),
        ArraySeqMajorBound([9, 8, 9, 8]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_B_1: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 2, 1]),
        ArraySeqMajorBound([3, 1, 3, 1]),
        ArraySeqMajorBound([4, 1, 4, 1]),
        ArraySeqMajorBound([5, 1, 5, 1]),
        ArraySeqMajorBound([5, 2, 5, 2]),
        ArraySeqMajorBound([6, 2, 6, 2]),
        ArraySeqMajorBound([6, 3, 6, 3]),
        ArraySeqMajorBound([7, 3, 7, 3]),
        ArraySeqMajorBound([7, 4, 7, 4]),
        ArraySeqMajorBound([8, 4, 8, 4]),
        ArraySeqMajorBound([8, 5, 8, 5]),
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
    pub const BOUNDS_B_2: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([2, 5, 2, 5]),
        ArraySeqMajorBound([2, 6, 2, 6]),
        ArraySeqMajorBound([3, 6, 3, 6]),
        ArraySeqMajorBound([3, 7, 3, 7]),
        ArraySeqMajorBound([4, 7, 4, 7]),
        ArraySeqMajorBound([4, 8, 4, 8]),
        ArraySeqMajorBound([5, 8, 5, 8]),
        ArraySeqMajorBound([5, 9, 5, 9]),
        ArraySeqMajorBound([6, 9, 6, 9]),
        ArraySeqMajorBound([7, 9, 7, 9]),
        ArraySeqMajorBound([8, 9, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_B_MERGE: [ArraySeqMajorBound; 17] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 2, 1]),
        ArraySeqMajorBound([3, 1, 2, 2]),
        ArraySeqMajorBound([4, 1, 2, 3]),
        ArraySeqMajorBound([5, 1, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([6, 2, 2, 6]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([7, 3, 3, 7]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([8, 4, 4, 8]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 5, 9]),
        ArraySeqMajorBound([8, 7, 6, 9]),
        ArraySeqMajorBound([8, 8, 7, 9]),
        ArraySeqMajorBound([8, 9, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_C_1: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 1, 2]),
        ArraySeqMajorBound([3, 1, 1, 3]),
        ArraySeqMajorBound([3, 2, 2, 3]),
        ArraySeqMajorBound([4, 2, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([5, 3, 3, 5]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([8, 4, 4, 8]),
        ArraySeqMajorBound([9, 4, 4, 9]),
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
    pub const BOUNDS_C_2: [ArraySeqMajorBound; 12] = [
        ArraySeqMajorBound([6, 1, 1, 6]),
        ArraySeqMajorBound([6, 2, 2, 6]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([7, 5, 5, 7]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([8, 7, 7, 8]),
        ArraySeqMajorBound([9, 7, 7, 9]),
        ArraySeqMajorBound([9, 8, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_C_MERGE: [ArraySeqMajorBound; 17] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 1, 2]),
        ArraySeqMajorBound([3, 1, 1, 3]),
        ArraySeqMajorBound([4, 1, 1, 4]),
        ArraySeqMajorBound([5, 1, 1, 5]),
        ArraySeqMajorBound([6, 1, 1, 6]),
        ArraySeqMajorBound([6, 2, 2, 6]),
        ArraySeqMajorBound([6, 3, 3, 6]),
        ArraySeqMajorBound([6, 4, 4, 6]),
        ArraySeqMajorBound([7, 4, 4, 7]),
        ArraySeqMajorBound([8, 4, 4, 8]),
        ArraySeqMajorBound([9, 4, 4, 9]),
        ArraySeqMajorBound([9, 5, 5, 9]),
        ArraySeqMajorBound([9, 6, 6, 9]),
        ArraySeqMajorBound([9, 7, 7, 9]),
        ArraySeqMajorBound([9, 8, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_D_1: [ArraySeqMajorBound; 5] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 1, 2]),
        ArraySeqMajorBound([3, 1, 1, 3]),
        ArraySeqMajorBound([3, 2, 2, 3]),
        ArraySeqMajorBound([4, 2, 2, 4]),
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
    pub const BOUNDS_D_2: [ArraySeqMajorBound; 5] = [
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([8, 7, 7, 8]),
        ArraySeqMajorBound([9, 7, 7, 9]),
        ArraySeqMajorBound([9, 8, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
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
    pub const BOUNDS_D_MERGE: [ArraySeqMajorBound; 17] = [
        ArraySeqMajorBound([1, 1, 1, 1]),
        ArraySeqMajorBound([2, 1, 1, 2]),
        ArraySeqMajorBound([3, 1, 1, 3]),
        ArraySeqMajorBound([3, 2, 2, 3]),
        ArraySeqMajorBound([4, 2, 2, 4]),
        ArraySeqMajorBound([5, 2, 2, 5]),
        ArraySeqMajorBound([6, 2, 2, 6]),
        ArraySeqMajorBound([7, 2, 2, 7]),
        ArraySeqMajorBound([8, 2, 2, 8]),
        ArraySeqMajorBound([8, 3, 3, 8]),
        ArraySeqMajorBound([8, 4, 4, 8]),
        ArraySeqMajorBound([8, 5, 5, 8]),
        ArraySeqMajorBound([8, 6, 6, 8]),
        ArraySeqMajorBound([8, 7, 7, 8]),
        ArraySeqMajorBound([9, 7, 7, 9]),
        ArraySeqMajorBound([9, 8, 8, 9]),
        ArraySeqMajorBound([9, 9, 9, 9]),
    ];
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_consts::*;

    use assert2::assert;

    #[test]
    pub fn test_advance_forward() {
        let mut cloud = Cloud::new(10, 10);
        cloud.set(2, 1, 1, 1, 1);
        cloud.ad_start = 2;
        cloud.ad_end = 2;

        (0..18).for_each(|_| cloud.advance_forward());

        let mut target_cloud = Cloud::new(10, 10);
        target_cloud.fill_rectangle(1, 1, 10, 10);

        assert!(cloud == target_cloud);
    }

    #[test]
    pub fn test_advance_reverse() {
        let mut cloud = Cloud::new(10, 10);
        cloud.set(20, 10, 10, 10, 10);
        cloud.ad_start = 20;
        cloud.ad_end = 20;

        (0..18).for_each(|_| cloud.advance_reverse());

        let mut target_cloud = Cloud::new(10, 10);
        target_cloud.fill_rectangle(1, 1, 10, 10);

        assert!(cloud == target_cloud);
    }

    #[test]
    pub fn test_fill_square() -> anyhow::Result<()> {
        let mut cloud = Cloud::new(10, 10);
        cloud.fill_rectangle(3, 3, 6, 6);
        let mut target_cloud = Cloud::new(10, 10);
        target_cloud.fill(&BOUNDS_SQUARE_3_3_6_6)?;
        assert!(cloud == target_cloud);

        Ok(())
    }

    #[test]
    pub fn test_trim_wings() -> anyhow::Result<()> {
        let mut cloud = Cloud::new(10, 10);
        cloud.fill(&BOUNDS_WINGS)?;

        let mut target_cloud = Cloud::new(10, 10);
        target_cloud.fill(&BOUNDS_WINGS_TRIMMED)?;

        cloud.trim_wings()?;
        assert!(cloud == target_cloud);

        // re-run it to make sure it does
        // nothing to already-trimmed cloud
        cloud.trim_wings()?;
        assert!(cloud == target_cloud);

        Ok(())
    }

    #[test]
    fn test_square_corners() -> anyhow::Result<()> {
        let mut cloud = Cloud::new(10, 10);
        cloud.fill(&BOUNDS_UNSQUARED)?;

        let mut target_cloud = Cloud::new(10, 10);
        target_cloud.fill(&BOUNDS_SQUARED)?;

        target_cloud.ad_start = 6;
        target_cloud.ad_end = 14;

        cloud.square_corners();
        assert!(cloud == target_cloud);

        // re-run it to make sure it does
        // nothing to already-squared cloud
        cloud.square_corners();
        assert!(cloud == target_cloud);

        Ok(())
    }

    #[test]
    fn test_merge() -> anyhow::Result<()> {
        let mut b1 = Cloud::new(10, 10);
        let mut b2 = Cloud::new(10, 10);
        let mut target_cloud = Cloud::new(10, 10);

        //--------------------------------
        b1.fill(&BOUNDS_A_1)?;
        b2.fill(&BOUNDS_A_2)?;
        target_cloud.fill(&BOUNDS_A_MERGE)?;

        b1.merge(&b2);
        assert!(b1 == target_cloud);

        // reset and invert merge order
        b1.fill(&BOUNDS_A_1)?;

        b2.merge(&b1);
        assert!(b2 == target_cloud);

        //--------------------------------
        b1.fill(&BOUNDS_B_1)?;
        b2.fill(&BOUNDS_B_2)?;
        target_cloud.fill(&BOUNDS_B_MERGE)?;

        b1.merge(&b2);
        assert!(b1 == target_cloud);

        // reset and invert merge order
        b1.fill(&BOUNDS_B_1)?;

        b2.merge(&b1);
        assert!(b2 == target_cloud);

        //--------------------------------
        b1.fill(&BOUNDS_C_1)?;
        b2.fill(&BOUNDS_C_2)?;
        target_cloud.fill(&BOUNDS_C_MERGE)?;

        b1.merge(&b2);
        assert!(b1 == target_cloud);

        // reset and invert merge order
        b1.fill(&BOUNDS_C_1)?;

        b2.merge(&b1);
        assert!(b2 == target_cloud);

        //--------------------------------
        b1.fill(&BOUNDS_D_1)?;
        b2.fill(&BOUNDS_D_2)?;
        target_cloud.fill(&BOUNDS_D_MERGE)?;

        b1.merge(&b2);

        assert!(b1 == target_cloud);

        // reset and invert merge order
        b1.fill(&BOUNDS_D_1)?;

        b2.merge(&b1);
        assert!(b2 == target_cloud);
        Ok(())
    }

    #[test]
    fn test_anti_diagonal_grow_up() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_target(6);
        assert!(a.right_target_idx == 5);
        assert!(a.right_profile_idx == 5);

        // make sure we grow properly
        a.grow_target(1);
        assert!(a.right_target_idx == 1);
        assert!(a.right_profile_idx == 9);

        // make sure we can't grow past 1
        a.grow_target(0);
        assert!(a.right_target_idx == 1);
        assert!(a.right_profile_idx == 9);
    }

    #[test]
    fn test_anti_diagonal_grow_down() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_down(4);
        assert!(a.left_target_idx == 5);
        assert!(a.left_profile_idx == 5);

        // make sure we grow properly
        a.grow_down(9);
        assert!(a.left_target_idx == 9);
        assert!(a.left_profile_idx == 1);
    }

    #[test]
    fn test_anti_diagonal_grow_left() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_left(6);
        assert!(a.left_target_idx == 5);
        assert!(a.left_profile_idx == 5);

        // make sure we grow properly
        a.grow_left(1);
        assert!(a.left_target_idx == 9);
        assert!(a.left_profile_idx == 1);

        // make sure we can't grow past 1
        a.grow_left(0);
        assert!(a.left_target_idx == 9);
        assert!(a.left_profile_idx == 1);
    }

    #[test]
    fn test_anti_diagonal_grow_right() {
        let mut a = AntiDiagonal::new(5, 5, 5, 5);

        // make sure we can't grow
        // in the wrong direction
        a.grow_right(4);
        assert!(a.right_target_idx == 5);
        assert!(a.right_profile_idx == 5);

        // make sure we grow properly
        a.grow_right(9);
        assert!(a.right_target_idx == 1);
        assert!(a.right_profile_idx == 9);
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
        let mut a = Cloud::new(10, 10);
        a.fill_rectangle(1, 1, 4, 4);

        let mut b = Cloud::new(10, 10);
        b.fill_rectangle(6, 6, 10, 10);

        // cloud should always intersect with themselves
        assert!(
            a.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 2, end: 8 })
        );
        assert!(
            b.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 12, end: 20 })
        );

        assert!(a.cloud_relationship(&b) == Relationship::Disjoint(Interval { start: 9, end: 11 }));
        assert!(b.cloud_relationship(&a) == Relationship::Disjoint(Interval { start: 9, end: 11 }));

        // --------------------------------------------------------
        // no overlap == touching
        let mut a = Cloud::new(10, 10);
        a.fill_rectangle(1, 1, 5, 6);

        let mut b = Cloud::new(10, 10);
        b.fill_rectangle(6, 5, 10, 10);

        // cloud should always intersect with themselves
        assert!(
            a.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 2, end: 11 })
        );
        assert!(
            b.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 11, end: 20 })
        );

        assert!(a.cloud_relationship(&b) == Relationship::Disjoint(Interval { start: 0, end: 0 }),);
        assert!(b.cloud_relationship(&a) == Relationship::Disjoint(Interval { start: 0, end: 0 }),);

        // --------------------------------------------------------
        // 1 cell overlap
        let mut a = Cloud::new(10, 10);
        a.fill_rectangle(1, 1, 5, 5);

        let mut b = Cloud::new(10, 10);
        b.fill_rectangle(5, 5, 10, 10);

        // cloud should always intersect with themselves
        assert!(
            a.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 2, end: 10 }),
        );
        assert!(
            b.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 10, end: 20 }),
        );

        assert!(
            a.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 10, end: 10 }),
        );
        assert!(
            b.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 10, end: 10 }),
        );

        // --------------------------------------------------------
        // many cell overlap
        let mut a = Cloud::new(10, 10);
        a.fill_rectangle(1, 1, 8, 8);

        let mut b = Cloud::new(10, 10);
        b.fill_rectangle(2, 2, 10, 10);

        // cloud should always intersect with themselves
        assert!(
            a.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 2, end: 16 }),
        );
        assert!(
            b.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 4, end: 20 }),
        );

        assert!(
            a.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 4, end: 16 }),
        );
        assert!(
            b.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 4, end: 16 }),
        );

        // --------------------------------------------------------
        // anti-diagonal overlap, not touching
        let mut a = Cloud::new(10, 10);
        a.fill_rectangle(1, 1, 8, 4);

        let mut b = Cloud::new(10, 10);
        b.fill_rectangle(1, 5, 8, 9);

        // cloud should always intersect with themselves
        assert!(
            a.cloud_relationship(&a) == Relationship::Intersecting(Interval { start: 2, end: 12 }),
        );
        assert!(
            b.cloud_relationship(&b) == Relationship::Intersecting(Interval { start: 6, end: 17 }),
        );

        assert!(a.cloud_relationship(&b) == Relationship::Disjoint(Interval { start: 0, end: 0 }),);
        assert!(b.cloud_relationship(&a) == Relationship::Disjoint(Interval { start: 0, end: 0 }),);
    }
}
