use crate::align::bounded::structs::bound::CloudBoundGroup;

#[derive(Default)]
pub struct RowBoundParams {
    pub target_start: usize,
    pub target_end: usize,
    pub row_capacity: usize,
    pub left_row_bounds: Vec<usize>,
    pub right_row_bounds: Vec<usize>,
}

impl RowBoundParams {
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
            self.left_row_bounds[bound.left_target_idx] =
                self.left_row_bounds[bound.left_target_idx].min(bound.left_profile_idx);

            self.right_row_bounds[bound.right_target_idx] =
                self.right_row_bounds[bound.right_target_idx].max(bound.right_profile_idx);
        }
    }
}
