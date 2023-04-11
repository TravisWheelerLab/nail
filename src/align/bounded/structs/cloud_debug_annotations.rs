use crate::align::bounded::structs::{CloudBoundGroup, RowBounds};
use crate::viz::SodaAnnotation;
use serde::Serialize;

#[derive(Serialize, Default)]
pub struct CloudDebugAnnotations {
    pub forward_bounds: (Vec<SodaAnnotation>, Vec<SodaAnnotation>),
    pub backward_bounds: (Vec<SodaAnnotation>, Vec<SodaAnnotation>),
    pub joined_bounds: (Vec<SodaAnnotation>, Vec<SodaAnnotation>),
    pub row_bounds: Vec<SodaAnnotation>,
}

impl CloudDebugAnnotations {
    pub fn add_forward_bounds(&mut self, bounds: &CloudBoundGroup) {
        Self::add_anti_diagonal_bounds(&mut self.forward_bounds, bounds);
    }

    pub fn add_backward_bounds(&mut self, bounds: &CloudBoundGroup) {
        Self::add_anti_diagonal_bounds(&mut self.backward_bounds, bounds);
    }

    pub fn add_joined_bounds(&mut self, bounds: &CloudBoundGroup) {
        Self::add_anti_diagonal_bounds(&mut self.joined_bounds, bounds);
    }

    pub fn add_row_bounds(&mut self, bounds: &RowBounds) {
        self.row_bounds = bounds.left_row_bounds[bounds.target_start..=bounds.target_end]
            .iter()
            .zip(&bounds.right_row_bounds[bounds.target_start..=bounds.target_end])
            .enumerate()
            .map(|(r, (s, e))| SodaAnnotation {
                id: format!("row-{}", r),
                start: *s,
                end: *e,
                row: r,
            })
            .collect();
    }

    fn add_anti_diagonal_bounds(
        bounds_ref: &mut (Vec<SodaAnnotation>, Vec<SodaAnnotation>),
        bounds: &CloudBoundGroup,
    ) {
        bounds_ref.0 = bounds.bounds[bounds.min_anti_diagonal_idx..=bounds.max_anti_diagonal_idx]
            .iter()
            .map(|b| SodaAnnotation {
                id: format!("left-{}", b.left_target_idx + b.left_profile_idx),
                start: b.left_profile_idx,
                end: b.left_profile_idx + 1,
                row: b.left_target_idx,
            })
            .collect();

        bounds_ref.1 = bounds.bounds[bounds.min_anti_diagonal_idx..=bounds.max_anti_diagonal_idx]
            .iter()
            .map(|b| SodaAnnotation {
                id: format!("right-{}", b.right_target_idx + b.right_profile_idx),
                start: b.right_profile_idx,
                end: b.right_profile_idx + 1,
                row: b.right_target_idx,
            })
            .collect();
    }
}
