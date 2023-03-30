use crate::align::bounded::structs::RowBoundParams;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::Profile;

#[derive(Default)]
pub struct DpMatrixSparse {
    pub target_length: usize,
    pub profile_length: usize,
    /// These point to the start of each "row block" in the `core_data` vector. They are absolute
    /// in the sense that they are not relative to the logical non-sparse matrix coordinate space.
    pub block_offsets: Vec<usize>,
    /// These indicate when a row "starts" in the logical non-sparse matrix coordinate space.
    pub row_offsets: Vec<usize>,
    pub core_data: Vec<f32>,
    pub special_data: Vec<f32>,
}

impl DpMatrixSparse {
    pub fn new(target_length: usize, profile_length: usize, row_bounds: &RowBoundParams) -> Self {
        let mut core_length = 0;
        // +1 for 0 row, + 1 for target_length + 1 row (serves as an end pointer)
        let mut block_offsets = vec![0; target_length + 2];
        // +1 for 0 row
        let mut row_offsets = vec![0; target_length + 1];
        for target_idx in row_bounds.target_start..=row_bounds.target_end {
            let row_length = (row_bounds.right_row_bounds[target_idx]
                - row_bounds.left_row_bounds[target_idx]
                // +1 for pad, +1 for subtraction across the interval
                + 2)
                // *3 for match, insert, delete cells
                * 3;

            core_length += row_length;
            block_offsets[target_idx + 1] = block_offsets[target_idx] + row_length;
            row_offsets[target_idx] = row_bounds.target_start;
        }

        let special_length = 5 * (target_length + 1);
        DpMatrixSparse {
            target_length,
            profile_length,
            block_offsets,
            row_offsets,
            core_data: vec![-f32::INFINITY; core_length],
            special_data: vec![-f32::INFINITY; special_length],
        }
    }

    pub fn resize(&mut self, new_target_length: usize, new_profile_length: usize) {
        todo!()
    }

    pub fn reset(&mut self) {
        todo!()
    }

    pub fn reuse(&mut self, new_target_length: usize, new_profile_length: usize) {
        todo!()
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
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset + 1]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset + 1] = value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset + 2]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        let block_offset = self.block_offsets[target_idx];
        let row_offset = self.row_offsets[target_idx];
        self.core_data[block_offset + (3 * profile_idx) - row_offset + 2] = value;
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
