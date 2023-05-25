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
pub struct CloudBoundGroup {
    pub bounds: Vec<CloudBound>,
    pub target_length: usize,
    pub profile_length: usize,
    pub size: usize,
    pub min_anti_diagonal_idx: usize,
    pub max_anti_diagonal_idx: usize,
}

impl CloudBoundGroup {
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
        for idx in self.min_anti_diagonal_idx..=self.max_anti_diagonal_idx {
            let bound = &self.bounds[idx];
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
    ///
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

    pub fn join_bounds(forward_bounds: &mut CloudBoundGroup, backward_bounds: &CloudBoundGroup) {
        // first check if the forward & backward bounds happened not to intersect
        if forward_bounds.max_anti_diagonal_idx < backward_bounds.min_anti_diagonal_idx {
            // TODO: **this is going to need to be rewritten**
            //       I think a better idea is going to be to interpolate two lines:
            //         1. one between the highest point on the forward & backward bounds
            //         2. one between the lowest point on the forward & backward bounds
            //
            // if they do not intersect, we need to interpolate
            //
            // we are going to do that by:
            //   1. selecting a "central point" on both the last
            //      forward bound and the first backward bound
            //   2. solving for the linear equation that defines
            //      the line between those two points
            //   3. computing the average anti-diagonal length
            //      across the anti-diagonals in both bound groups
            //   4. fill in between the bound groups with anti-diagonals
            //      of that average length that are roughly centered
            //      on the line
            //
            let last_forward = forward_bounds.get_last();
            let forward_x = ((last_forward.right_profile_idx as f32
                + last_forward.left_profile_idx as f32)
                / 2.0)
                .floor();
            let forward_y = last_forward.anti_diagonal_idx() as f32 - forward_x;

            let first_backward = backward_bounds.get_first();
            let backward_x = ((first_backward.right_profile_idx as f32
                + first_backward.left_profile_idx as f32)
                / 2.0)
                .ceil();
            let backward_y = first_backward.anti_diagonal_idx() as f32 - backward_x;

            let slope = (backward_y - forward_y) / (backward_x - forward_x);
            let intercept = forward_y - slope * forward_x;
            let line_equation = |x: f32| slope * x + intercept;

            let cloud_size_sum = forward_bounds.cloud_size() + backward_bounds.cloud_size();
            let num_anti_diagonals =
                forward_bounds.num_anti_diagonals() + backward_bounds.num_anti_diagonals();
            // integer division is truncated, which is probably what we want
            let avg_anti_diagonal_length = cloud_size_sum / num_anti_diagonals;

            // I have elected to make this a closure since it's convenient
            // to have it capture the kine_equation closure, and I don't think
            // we'll ever call this from outside of this function
            let bound_fn = |center_profile_idx: usize| {
                let center_target_idx = line_equation(center_profile_idx as f32).round() as usize;
                let anti_diagonal_idx = center_profile_idx + center_target_idx;

                // note: if the avg_anti_diagonal_length is even, this
                //       will cause us to fill with +1 of that length
                let left_profile_idx = center_profile_idx - avg_anti_diagonal_length / 2;
                let right_profile_idx = center_profile_idx + avg_anti_diagonal_length / 2;
                let left_target_idx = anti_diagonal_idx - left_profile_idx;
                let right_target_idx = anti_diagonal_idx - right_profile_idx;
                CloudBound {
                    left_target_idx,
                    left_profile_idx,
                    right_target_idx,
                    right_profile_idx,
                }
            };

            let first_bound = bound_fn(forward_x as usize);
            forward_bounds.set(
                first_bound.anti_diagonal_idx() + 1,
                first_bound.left_target_idx,
                first_bound.left_profile_idx + 1,
                first_bound.right_target_idx + 1,
                first_bound.right_profile_idx,
            );

            // fill in the missing bounds using the line equation:
            //   - iterate across the profile indices that span the missing bounds
            //   - we'll think of the profile index as the x value in the line equation
            //   - using the profile index as the input to the line equation, we can get
            //     a corresponding target index (y value)
            //   - those coordinates (profile index, target index) will be the center cell
            //     along the current anti-diagonal that we are filling
            //   - then we can just extend that bound out from the center to produce an
            //     anti-diagonal equal to the average anti-diagonal length
            let profile_start = forward_x as usize + 1;
            let profile_end = backward_x as usize - 1;
            for center_profile_idx in profile_start..=profile_end {
                let bound = bound_fn(center_profile_idx);

                forward_bounds.set(
                    bound.anti_diagonal_idx(),
                    bound.left_target_idx,
                    bound.left_profile_idx,
                    bound.right_target_idx,
                    bound.right_profile_idx,
                );

                forward_bounds.set(
                    bound.anti_diagonal_idx() + 1,
                    bound.left_target_idx,
                    bound.left_profile_idx + 1,
                    bound.right_target_idx + 1,
                    bound.right_profile_idx,
                );
            }

            // finally just tack on the backward bounds
            for anti_diagonal_idx in
                backward_bounds.min_anti_diagonal_idx..=backward_bounds.max_anti_diagonal_idx
            {
                let bound = backward_bounds.get(anti_diagonal_idx);
                forward_bounds.set(
                    anti_diagonal_idx,
                    bound.left_target_idx,
                    bound.left_profile_idx,
                    bound.right_target_idx,
                    bound.right_profile_idx,
                )
            }
        } else {
            // if they do intersect, we can join them by taking
            // the longest anti-diagonal at each index
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
}
