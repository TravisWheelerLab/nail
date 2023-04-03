use crate::align::bounded::structs::RowBoundParams;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::Profile;

#[derive(Default)]
pub struct DpMatrixSparse {
    pub target_length: usize,
    pub profile_length: usize,
    pub target_start: usize,
    pub target_end: usize,
    /// These point to the start of each "row block" in the `core_data` vector. They are absolute
    /// in the sense that they are not relative to the logical non-sparse matrix coordinate space.
    pub block_offsets: Vec<usize>,
    /// These indicate when a row "starts" in the logical non-sparse matrix coordinate space.
    pub row_start_offsets: Vec<usize>,
    pub core_data: Vec<f32>,
    pub special_data: Vec<f32>,
}

impl DpMatrixSparse {
    pub fn new(target_length: usize, profile_length: usize, row_bounds: &RowBoundParams) -> Self {
        let mut matrix = DpMatrixSparse::default();
        matrix.reuse(target_length, profile_length, row_bounds);
        matrix
    }

    pub fn reset(&mut self) {
        for core_idx in 0..self.block_offsets[self.target_end + 1] {
            self.core_data[core_idx] = -f32::INFINITY;
        }

        for special_idx in 0..(5 * (self.target_length + 1)) {
            self.special_data[special_idx] = -f32::INFINITY;
        }
    }

    pub fn reuse(
        &mut self,
        new_target_length: usize,
        new_profile_length: usize,
        row_bounds: &RowBoundParams,
    ) {
        let mut core_length = 0;
        // +1 for 0 row, + 1 for target_length + 1 row (serves as an end pointer)
        let mut block_offsets = vec![0; new_target_length + 2];
        // +1 for 0 row
        let mut row_offsets = vec![0; new_target_length + 1];

        // iterate across the row bounds to compute
        // how many cells we need for the core data
        for target_idx in row_bounds.target_start..=row_bounds.target_end {
            let row_length = (row_bounds.right_row_bounds[target_idx]
                - row_bounds.left_row_bounds[target_idx]
                // +1 for pad, +1 for subtraction across the interval
                + 2)
                // *3 for match, insert, delete cells
                * 3;

            core_length += row_length;
            block_offsets[target_idx + 1] = block_offsets[target_idx] + row_length;
            row_offsets[target_idx] = row_bounds.left_row_bounds[target_idx] - 1;
        }

        let new_special_length = 5 * (new_target_length + 1);

        if core_length > self.core_data.len() {
            self.core_data.resize(core_length, -f32::INFINITY);
        }

        if new_special_length > self.special_data.len() {
            self.special_data.resize(new_special_length, -f32::INFINITY);
        }

        self.row_start_offsets = row_offsets;
        self.block_offsets = block_offsets;
        self.target_length = new_target_length;
        self.profile_length = new_profile_length;

        self.target_start = row_bounds.target_start;
        self.target_end = row_bounds.target_end;

        self.reset();
    }
}

impl DpMatrix for DpMatrixSparse {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    // TODO: experiment with index computing functions

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset)))]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset)))] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset))) + 1]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset))) + 1] = value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset))) + 2]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_start_offsets[target_idx];
        self.core_data[block_offset + (3 * (profile_idx.saturating_sub(row_offset))) + 2] = value;
    }

    #[inline]
    fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx]
    }

    #[inline]
    fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx] = value;
    }
}
