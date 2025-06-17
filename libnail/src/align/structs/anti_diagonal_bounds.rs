use anyhow::{anyhow, bail};

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::io::{stdout, Write};
use std::ops::{Index, IndexMut};

use crate::util::{MaxAssign, MinAssign, VecUtils};

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

#[derive(Default, Clone, Copy, PartialEq)]
pub struct Cell {
    pub prf_idx: usize,
    pub seq_idx: usize,
}

impl Cell {
    pub(crate) fn from_seq_major(src: &[usize]) -> Self {
        Self {
            prf_idx: src[1],
            seq_idx: src[0],
        }
    }

    pub(crate) fn from_profile_major(src: &[usize]) -> Self {
        Self {
            prf_idx: src[0],
            seq_idx: src[1],
        }
    }

    pub fn idx(&self) -> usize {
        self.seq_idx + self.prf_idx
    }

    pub fn m_cell(&self) -> CoreCell {
        (CoreState::M(self.prf_idx), self.seq_idx)
    }

    pub fn i_cell(&self) -> CoreCell {
        (CoreState::I(self.prf_idx), self.seq_idx)
    }

    pub fn d_cell(&self) -> CoreCell {
        (CoreState::D(self.prf_idx), self.seq_idx)
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
        write!(f, "(p{}, s{})", self.prf_idx, self.seq_idx)
    }
}

#[allow(dead_code)]
pub(crate) struct ArrayProfileMajorBound([usize; 4]);

#[allow(dead_code)]
pub(crate) struct ArraySeqMajorBound([usize; 4]);

#[derive(Copy, Clone, PartialEq)]
pub struct Bound(pub Cell, pub Cell);

impl Default for Bound {
    fn default() -> Self {
        Self(
            Cell {
                prf_idx: 0,
                seq_idx: 0,
            },
            Cell {
                prf_idx: 0,
                seq_idx: 0,
            },
        )
    }
}

impl Bound {
    pub fn idx(&self) -> usize {
        debug_assert!(self.0.idx() == self.1.idx());
        self.0.idx()
    }

    pub fn iter(&self) -> BoundIter {
        BoundIter {
            profile_start: self.0.prf_idx,
            profile_end: self.1.prf_idx,
            seq_start: self.0.seq_idx,
            seq_end: self.1.seq_idx,
        }
    }

    pub fn len(&self) -> usize {
        self.0.prf_idx - self.1.prf_idx + 1
    }

    pub fn is_empty(&self) -> bool {
        self.0.prf_idx < self.1.prf_idx
    }

    pub fn is_valid(&self) -> bool {
        self.0.prf_idx >= self.1.prf_idx
            && self.1.seq_idx >= self.0.seq_idx
            && (self.0.prf_idx + self.0.seq_idx) == (self.1.prf_idx + self.1.seq_idx)
    }

    pub fn intersects(&self, other: &Self) -> bool {
        if self.idx() != other.idx() {
            false
        } else {
            self.1.prf_idx <= other.0.prf_idx && self.0.prf_idx >= other.1.prf_idx
        }
    }

    pub fn is_default(&self) -> bool {
        *self == Self::default()
    }

    pub fn replace_with(&mut self, other: &Self) {
        self.0 = other.0;
        self.1 = other.1;
    }

    pub fn grow_down_seq(&mut self, seq_idx: usize) {
        let idx = self.idx();
        self.0.seq_idx.min_assign(seq_idx.max(1));
        self.0.prf_idx = idx - self.0.seq_idx;
    }

    pub fn grow_up_seq(&mut self, seq_idx: usize) {
        let idx = self.idx();
        self.1.seq_idx.max_assign(seq_idx);
        self.0.seq_idx = idx - self.1.seq_idx;
    }

    pub fn grow_up_prf(&mut self, prf_idx: usize) {
        let idx = self.idx();
        self.0.seq_idx.min_assign(prf_idx.max(1));
        self.1.seq_idx = idx - self.0.seq_idx;
    }

    pub fn grow_down_prf(&mut self, prf_idx: usize) {
        let idx = self.idx();
        self.0.prf_idx.max_assign(prf_idx);
        self.0.seq_idx = idx - self.0.prf_idx;
    }

    pub fn merge(&mut self, other: &Self) {
        self.0.prf_idx.max_assign(other.0.prf_idx);
        self.0.seq_idx.min_assign(other.0.seq_idx);
        self.1.prf_idx.min_assign(other.1.prf_idx);
        self.1.seq_idx.max_assign(other.1.seq_idx);
    }

    #[allow(dead_code)]
    pub(crate) fn from_seq_major(src: &[usize]) -> Self {
        Bound(
            Cell::from_seq_major(&src[2..4]),
            Cell::from_seq_major(&src[0..2]),
        )
    }

    #[allow(dead_code)]
    pub(crate) fn from_profile_major(src: &[usize]) -> Self {
        Bound(
            Cell::from_profile_major(&src[0..2]),
            Cell::from_profile_major(&src[2..4]),
        )
    }
}

impl Display for Bound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(p{}, s{}), (p{}, s{}) | AD: {}",
            self.0.prf_idx,
            self.0.seq_idx,
            self.1.prf_idx,
            self.1.seq_idx,
            self.idx()
        )
    }
}

impl Debug for Bound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(p{}, s{}), (p{}, s{}) | AD: {}",
            self.0.prf_idx,
            self.0.seq_idx,
            self.1.prf_idx,
            self.1.seq_idx,
            self.idx()
        )
    }
}

pub struct BoundIter {
    profile_start: usize,
    profile_end: usize,
    seq_start: usize,
    seq_end: usize,
}

impl Debug for BoundIter {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mid = if self.profile_end <= self.profile_start {
            "-".repeat(self.profile_start - self.profile_end + 1)
        } else {
            "".to_string()
        };
        writeln!(
            f,
            "p{}>{mid}<{}\ns{}>{mid}<{}",
            self.profile_start, self.profile_end, self.seq_start, self.seq_end
        )
    }
}

impl Iterator for BoundIter {
    type Item = Cell;

    fn next(&mut self) -> Option<Self::Item> {
        if self.seq_start <= self.seq_end {
            let cell = Cell {
                prf_idx: self.profile_start,
                seq_idx: self.seq_start,
            };
            // saturating sub catches edge case (p0, s<n>)
            self.profile_start = self.profile_start.saturating_sub(1);
            self.seq_start += 1;
            Some(cell)
        } else {
            None
        }
    }
}

impl DoubleEndedIterator for BoundIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.profile_start >= self.profile_end {
            let cell = Cell {
                prf_idx: self.profile_end,
                seq_idx: self.seq_end,
            };
            self.profile_end += 1;
            // saturating sub catches edge case (p<n>, s0)
            self.seq_end = self.seq_end.saturating_sub(1);
            Some(cell)
        } else {
            None
        }
    }
}

pub trait BoundInterpretable {
    fn interpret(&self) -> Bound;
}

impl BoundInterpretable for Bound {
    fn interpret(&self) -> Bound {
        *self
    }
}

impl BoundInterpretable for ArrayProfileMajorBound {
    fn interpret(&self) -> Bound {
        Bound::from_profile_major(&self.0)
    }
}

impl BoundInterpretable for ArraySeqMajorBound {
    fn interpret(&self) -> Bound {
        Bound::from_seq_major(&self.0)
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
    pub bounds: Vec<Bound>,
    pub seq_len: usize,
    pub prf_len: usize,
    pub size: usize,
    pub ad_start: usize,
    pub ad_end: usize,
}

impl Default for Cloud {
    fn default() -> Self {
        Self {
            bounds: vec![Bound::default()],
            seq_len: 0,
            prf_len: 0,
            size: 1,
            ad_start: 0,
            ad_end: 0,
        }
    }
}

pub struct Ad(pub usize);
pub struct AdOffset(pub usize);

impl Index<Ad> for Cloud {
    type Output = Bound;

    fn index(&self, index: Ad) -> &Self::Output {
        &self.bounds[index.0]
    }
}

impl IndexMut<Ad> for Cloud {
    fn index_mut(&mut self, index: Ad) -> &mut Self::Output {
        self.ad_start = self.ad_start.min(index.0);
        self.ad_end = self.ad_end.max(index.0);
        &mut self.bounds[index.0]
    }
}

impl Index<AdOffset> for Cloud {
    type Output = Bound;

    fn index(&self, index: AdOffset) -> &Self::Output {
        &self.bounds[self.ad_start + index.0]
    }
}

impl IndexMut<AdOffset> for Cloud {
    fn index_mut(&mut self, index: AdOffset) -> &mut Self::Output {
        &mut self.bounds[self.ad_start + index.0]
    }
}

impl Cloud {
    pub fn new(seq_len: usize, prf_len: usize) -> Self {
        let mut cloud = Self::default();
        cloud.reuse(seq_len, prf_len);
        cloud
    }

    pub fn iter(&self) -> std::slice::Iter<Bound> {
        self.bounds[self.ad_start..=self.ad_end].iter()
    }

    pub fn first(&self) -> &Bound {
        &self.bounds[self.ad_start]
    }

    pub fn last(&self) -> &Bound {
        &self.bounds[self.ad_end]
    }

    pub fn resize(&mut self, new_size: usize) {
        self.bounds.grow_or_shrink(new_size, Bound::default());
    }

    pub fn reset(&mut self) {
        self.bounds.reset(Bound::default());
    }

    pub fn reuse(&mut self, seq_len: usize, prf_len: usize) {
        let new_size = seq_len + prf_len + 2;

        self.bounds.resize_and_reset(new_size, Bound::default());

        self.ad_start = new_size;
        self.ad_end = 0;

        self.seq_len = seq_len;
        self.prf_len = prf_len;
    }

    pub fn append(&mut self, bound: Bound) {
        let idx = bound.idx();
        self[Ad(idx)] = bound;
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

        bounds
            .into_iter()
            .for_each(|bound| self[Ad(bound.idx())] = bound);

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
        self.seq_len
    }

    pub fn profile_len(&self) -> usize {
        self.prf_len
    }

    pub fn valid(&self) -> bool {
        self.bounds.iter().all(|bound| bound.is_valid())
    }

    pub fn validity(&self) -> CloudValidity {
        let mut out_of_range = self.bounds[0..self.ad_start]
            .iter()
            .chain(self.bounds[self.ad_end + 1..].iter());

        let mut in_range = self.bounds[self.ad_start..=self.ad_end].iter();

        if let Some(ad) = out_of_range.find(|b| !b.is_default()) {
            CloudValidity::BoundOutOfRange {
                ad_start: self.ad_start,
                ad_end: self.ad_end,
                bound: *ad,
            }
        } else if let Some(ad) = in_range.find(|b| !b.is_valid()) {
            CloudValidity::BoundInvalid { bound: *ad }
        } else {
            CloudValidity::Valid
        }
    }

    /// Get the total number of cells that exist within the cloud boundaries.
    pub fn cloud_size(&self) -> usize {
        let mut cloud_size = 0usize;
        for bound in self.bounds[self.ad_start..=self.ad_end].iter() {
            cloud_size += bound.len()
        }
        cloud_size
    }

    /// Get the number of anti-diagonals defined in the cloud.
    pub fn num_anti_diagonals(&self) -> usize {
        self.ad_end - self.ad_start + 1
    }

    /// This removes all of the protruding regions in the cloud that are
    /// unreachable from a traceback that traverses the entire cloud.
    pub fn trim_wings(&mut self) -> anyhow::Result<()> {
        for idx in self.ad_start + 1..=self.ad_end {
            let prev = self[Ad(idx - 1)];
            let bound = &mut self[Ad(idx)];

            // the number of cells the right bound
            // is above the previous right bound
            let right_distance = prev.0.seq_idx.saturating_sub(bound.0.seq_idx);

            // the number of cells the left bound
            // is below the previous left bound
            let left_distance = (prev.1.prf_idx).saturating_sub(bound.1.prf_idx);

            bound.1.seq_idx -= left_distance;
            bound.1.prf_idx += left_distance;
            bound.0.seq_idx += right_distance;
            bound.0.prf_idx -= right_distance;

            if !bound.is_valid() {
                bail!("forward cloud trimming failed: invalid bound: {bound:?}");
            } else if bound.idx() - 1 != prev.idx() {
                println!("{prev}");
                println!("{bound}");
                println!("{left_distance} | {right_distance}");
                panic!(
                    "forward cloud trimming failed: AD index mismatch: {} | {}",
                    bound.idx(),
                    prev.idx()
                );
            }
        }

        for idx in (self.ad_start..self.ad_end).rev() {
            let next = self[Ad(idx + 1)];
            let bound = &mut self[Ad(idx)];

            // the number of cells the right bound
            // is above the previous right bound
            let right_distance = bound.0.prf_idx.saturating_sub(next.0.prf_idx);

            // the number of cells the left bound
            // is below the previous left bound
            let left_distance = (bound.1.seq_idx).saturating_sub(next.1.seq_idx);

            bound.1.seq_idx -= left_distance;
            bound.1.prf_idx += left_distance;
            bound.0.seq_idx += right_distance;
            bound.0.prf_idx -= right_distance;

            if !bound.is_valid() {
                bail!("reverse cloud trimming failed: invalid bound: {bound:?}");
            } else if bound.idx() + 1 != next.idx() {
                panic!(
                    "reverse cloud trimming failed: AD index mismatch: {} | {}",
                    bound.idx(),
                    next.idx()
                );
            }
        }

        Ok(())
    }

    pub fn square_corners(&mut self) {
        let first = *self.first();
        let distance = first.1.seq_idx - first.0.seq_idx;
        let ad_start = self.ad_start;
        (1..=distance)
            .map(|offset| (offset, ad_start - offset))
            .for_each(|(offset, idx)| {
                self.bounds[idx] = Bound(
                    Cell {
                        prf_idx: first.0.prf_idx - offset,
                        seq_idx: first.0.seq_idx,
                    },
                    Cell {
                        prf_idx: first.1.prf_idx,
                        seq_idx: first.1.seq_idx - offset,
                    },
                );
            });

        let last = *self.last();
        let distance = last.1.seq_idx - last.0.seq_idx;
        let ad_end = self.ad_end;
        (1..=distance)
            .map(|offset| (offset, ad_end + offset))
            .for_each(|(offset, idx)| {
                self[Ad(idx)] = Bound(
                    Cell {
                        prf_idx: last.0.prf_idx,
                        seq_idx: last.0.seq_idx + offset,
                    },
                    Cell {
                        prf_idx: last.1.prf_idx + offset,
                        seq_idx: last.1.seq_idx,
                    },
                );
            });
    }

    pub fn advance_forward(&mut self) {
        let last = self.bounds[self.ad_end];

        // this is going to be 1 if the anti-diagonal is going
        // to try to move past the target length, and 0 otherwise
        let prf_add = ((last.1.seq_idx + 1) > self.seq_len) as usize;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move past the profile length, and 0 otherwise
        let seq_add = ((last.0.prf_idx + 1) > self.prf_len) as usize;

        self.bounds[self.ad_end + 1].0.prf_idx = (last.0.prf_idx + 1).min(self.prf_len);
        self.bounds[self.ad_end + 1].0.seq_idx = last.0.seq_idx + seq_add;
        self.bounds[self.ad_end + 1].1.prf_idx = last.1.prf_idx + prf_add;
        self.bounds[self.ad_end + 1].1.seq_idx = (last.1.seq_idx + 1).min(self.seq_len);

        self.ad_end += 1;
    }

    pub fn advance_reverse(&mut self) {
        let first = self.bounds[self.ad_start];

        // this is going to be 1 if the anti-diagonal is going
        // to try to move before target index 1, and 0 otherwise
        let prf_sub = ((first.0.seq_idx.saturating_sub(1)) < 1) as usize;

        // this is going to be 1 if the anti-diagonal is going
        // to try to move before profle index 1, and 0 otherwise
        let seq_sub = ((first.1.prf_idx.saturating_sub(1)) < 1) as usize;

        self.bounds[self.ad_start - 1].0.prf_idx = first.0.prf_idx - prf_sub;
        self.bounds[self.ad_start - 1].0.seq_idx = (first.0.seq_idx - 1).max(1);
        self.bounds[self.ad_start - 1].1.prf_idx = (first.1.prf_idx - 1).max(1);
        self.bounds[self.ad_start - 1].1.seq_idx = first.1.seq_idx - seq_sub;

        self.ad_start -= 1;
    }

    pub fn fill_rectangle(
        &mut self,
        seq_start: usize,
        prf_start: usize,
        seq_end: usize,
        prf_end: usize,
    ) {
        self.reset();

        let seq_distance = seq_end - seq_start;
        let prf_distance = prf_end - prf_start;

        let ad_start = seq_start + prf_start;
        let ad_end = seq_end + prf_end;

        for idx in ad_start..=ad_end {
            // relative_idx is the number of antidiagonals we have moved
            // forward relative to the starting antidiagonal index
            let relative_idx = idx - ad_start;

            self.bounds[idx].0.seq_idx = seq_start + relative_idx.saturating_sub(prf_distance);
            self.bounds[idx].0.prf_idx = prf_end.min(prf_start + relative_idx);

            self.bounds[idx].1.seq_idx = seq_end.min(seq_start + relative_idx);
            self.bounds[idx].1.prf_idx = prf_start + relative_idx.saturating_sub(seq_distance);
        }

        self.ad_start = ad_start;
        self.ad_end = ad_end;
    }

    pub fn bounding_box(&self) -> BoundingBox {
        let mut bbox = BoundingBox::default();

        self.bounds.iter().for_each(|bound| {
            bbox.target_start = bbox.target_start.min(bound.1.seq_idx);
            bbox.target_end = bbox.target_start.max(bound.0.seq_idx);
            bbox.profile_start = bbox.profile_start.min(bound.1.prf_idx);
            bbox.profile_end = bbox.profile_start.max(bound.0.prf_idx);
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

            self_slice.iter_mut().zip(other_slice).for_each(|(b1, b2)| {
                b1.replace_with(b2);
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
                // merge along the intersecting anti-diagonals
                self.bounds[interval.start..=interval.end]
                    .iter_mut()
                    .zip(other.bounds[interval.start..=interval.end].iter())
                    .for_each(|(b1, b2)| b1.merge(b2));
            }
            Relationship::Disjoint(_) => panic!("tried to merge disjoint cloud"),
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
                for idx in interval.start..=interval.end {
                    let self_bound = &self[Ad(idx)];
                    let other_bound = &other[Ad(idx)];

                    if self_bound.intersects(other_bound) {
                        maybe_start = Some(idx);
                        break;
                    }
                }

                match maybe_start {
                    // if we've found the start of a proper intersection,
                    // we'll check to find the end of the intersection
                    Some(start) => {
                        let mut end = start;

                        for idx in (interval.start..=interval.end).rev() {
                            let self_bound = &self[Ad(idx)];
                            let other_bound = &other[Ad(idx)];

                            if self_bound.intersects(other_bound) {
                                end = idx;
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
        let mut img = vec![vec![0; self.seq_len + 2]; self.prf_len + 2];

        self.iter().for_each(|b| {
            b.iter().for_each(|c| img[c.prf_idx][c.seq_idx] += 1);
        });
        if let Some(other) = other {
            if self.prf_len != other.prf_len || self.seq_len != other.seq_len {
                bail!("dimension mismatch")
            }

            other.iter().for_each(|b| {
                b.iter().for_each(|c| img[c.prf_idx][c.seq_idx] += 2);
            });
        };
        Ok(img)
    }

    #[allow(dead_code)]
    pub fn image(&self) -> CloudImage {
        let img_bytes = self.vec_image(None).unwrap();
        CloudImage {
            img_bytes,
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

pub struct CloudImage {
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
            let path = path.as_ref();
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
                        bail!("cloud cell error: {val}");
                    }
                }
            }

            img.save(path)
                .with_context(|| format!("failed to write image to: {}", path.to_string_lossy()))
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
    pub fn test_bound_iter() -> anyhow::Result<()> {
        let bounds: Vec<Bound> = BOUNDS_A_MERGE.iter().map(|b| b.interpret()).collect();

        // note: sentinel [0,0] added to each vec to make
        // sure Bound::iter() stops where it's supposd to
        let cells_by_bound: Vec<Vec<Cell>> = [
            vec![[0, 0], [1, 1], [0, 0]],
            vec![[0, 0], [2, 1], [0, 0]],
            vec![[0, 0], [3, 1], [2, 2], [0, 0]],
            vec![[0, 0], [4, 1], [3, 2], [2, 3], [0, 0]],
            vec![[0, 0], [4, 2], [3, 3], [2, 4], [0, 0]],
            vec![[0, 0], [5, 2], [4, 3], [3, 4], [2, 5], [0, 0]],
            vec![[0, 0], [5, 3], [4, 4], [3, 5], [0, 0]],
            vec![[0, 0], [6, 3], [5, 4], [4, 5], [3, 6], [0, 0]],
            vec![[0, 0], [6, 4], [5, 5], [4, 6], [0, 0]],
            vec![[0, 0], [7, 4], [6, 5], [5, 6], [4, 7], [0, 0]],
            vec![[0, 0], [7, 5], [6, 6], [5, 7], [0, 0]],
            vec![[0, 0], [8, 5], [7, 6], [6, 7], [5, 8], [0, 0]],
            vec![[0, 0], [8, 6], [7, 7], [6, 8], [0, 0]],
            vec![[0, 0], [8, 7], [7, 8], [6, 9], [0, 0]],
            vec![[0, 0], [8, 8], [7, 9], [0, 0]],
            vec![[0, 0], [8, 9], [0, 0]],
            vec![[0, 0], [9, 9], [0, 0]],
        ]
        .iter()
        .map(|x| x.iter().map(|y| Cell::from_profile_major(y)).collect())
        .collect();

        bounds
            .iter()
            .zip(cells_by_bound)
            .for_each(|(bound, cells)| {
                assert!(bound.iter().count() == bound.iter().rev().count());

                bound
                    .iter()
                    // skip the leading sentinel (0, 0)
                    .zip(cells.iter().skip(1))
                    .for_each(|(c1, c2)| assert!(c1 == *c2));

                bound
                    .iter()
                    .rev()
                    // skip the trailing sentinel (0, 0)
                    .zip(cells.iter().rev().skip(1))
                    .for_each(|(c1, c2)| assert!(c1 == *c2));
            });

        Ok(())
    }

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
