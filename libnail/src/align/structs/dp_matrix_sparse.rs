use crate::align::structs::RowBounds;
use crate::structs::Profile;
use std::cell::RefMut;

use super::DpMatrix;

#[derive(Default, Clone)]
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
    pub fn new(target_length: usize, profile_length: usize, row_bounds: &RowBounds) -> Self {
        let mut matrix = DpMatrixSparse::default();
        matrix.reuse(target_length, profile_length, row_bounds);
        matrix
    }

    pub fn reset(&mut self) {
        // TODO: this resets more memory than is necessary
        for value in self.core_data.iter_mut() {
            *value = -f32::INFINITY;
        }

        for value in self.special_data.iter_mut() {
            *value = -f32::INFINITY;
        }
    }

    pub fn reuse(
        &mut self,
        new_target_length: usize,
        new_profile_length: usize,
        row_bounds: &RowBounds,
    ) {
        // we end up needing an extra 3 slots to make everything fit
        let mut core_length = 3;
        // +1 for the zero row, + 1 for the (target_length + 1) row
        let mut block_offsets = vec![0; new_target_length + 2];
        let mut row_offsets = vec![0; new_target_length + 2];

        // the first row in the sparse matrix is a row consisting entirely of pad cells
        // it should cover dependencies to the (target_start) row
        let first_row_idx = row_bounds.target_start - 1;
        let second_row_idx = row_bounds.target_start;

        // -1 to shift to the left one
        let first_row_start = row_bounds.left_row_bounds[second_row_idx] - 1;
        let first_row_end = row_bounds.right_row_bounds[second_row_idx] - 1;
        // +1 for subtraction across the interval
        let first_row_length = first_row_end - first_row_start + 1;

        row_offsets[first_row_idx] = first_row_start;
        // *3 for match, insert, delete cells
        core_length += first_row_length * 3;

        // the first row is always going to have a block
        // offset of 0, since it's the first thing stored
        block_offsets[first_row_idx] = 0;
        block_offsets[second_row_idx] = first_row_length * 3;

        for row_idx in row_bounds.target_start..=row_bounds.target_end {
            // -1 to add a pad cell
            let row_start = row_bounds.left_row_bounds[row_idx] - 1;
            let row_end = row_bounds.right_row_bounds[row_idx];

            // +1 for subtraction across the interval
            let row_length = row_end - row_start + 1;

            row_offsets[row_idx] = row_start;
            // *3 for match, insert, delete cells
            core_length += row_length * 3;

            // we set the block offset for the next row to avoid
            // having to keep track of the individual block lengths
            block_offsets[row_idx + 1] = block_offsets[row_idx] + (row_length * 3);
        }

        let last_row_idx = row_bounds.target_end + 1;
        let second_to_last_row_idx = row_bounds.target_end;

        let last_row_start = row_bounds.left_row_bounds[second_to_last_row_idx];
        // +1 to shift to the right one
        let last_row_end = row_bounds.right_row_bounds[second_to_last_row_idx] + 1;

        // +1 for subtraction across the interval
        let last_row_length = last_row_end - last_row_start + 1;

        row_offsets[last_row_idx] = last_row_start;
        // *3 for match, insert, delete cells
        core_length += last_row_length * 3;

        // the last block offset (target_end + 1) points to
        // the end of the last valid row in the matrix
        let last_block_offset_length =
            block_offsets[row_bounds.target_end + 1] + last_row_length * 3;

        let last_block_offset_idx = (row_bounds.target_end + 1).min(block_offsets.len());

        // we want to set everything past (target_end + 1) to point to the last block offset
        // this guarantees we can index into the entire logical matrix coordinate space
        for block_offset in
            block_offsets[(row_bounds.target_end + 1)..last_block_offset_idx].iter_mut()
        {
            *block_offset = last_block_offset_length;
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

    #[inline]
    fn match_idx(&self, target_idx: usize, profile_idx: usize) -> usize {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);

        let block_offset = self.block_offsets[target_idx];
        let next_block_offset = self.block_offsets[target_idx + 1];
        let row_offset = self.row_start_offsets[target_idx];

        // the saturating sub prevents moving left out of the block
        // the min with the next block offset prevents moving right out of the block
        (block_offset + (3 * (profile_idx.saturating_sub(row_offset)))).min(next_block_offset)
    }

    #[inline]
    fn insert_idx(&self, target_idx: usize, profile_idx: usize) -> usize {
        self.match_idx(target_idx, profile_idx) + 1
    }

    #[inline]
    fn delete_idx(&self, target_idx: usize, profile_idx: usize) -> usize {
        self.match_idx(target_idx, profile_idx) + 2
    }
}

impl DpMatrix for DpMatrixSparse {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.match_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.match_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.insert_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.insert_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.delete_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.delete_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
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

// TODO: I think we can get rid of this impl now
impl DpMatrix for RefMut<'_, DpMatrixSparse> {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.match_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.match_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.insert_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.insert_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        let idx = self.delete_idx(target_idx, profile_idx);
        self.core_data[idx]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        let idx = self.delete_idx(target_idx, profile_idx);
        self.core_data[idx] = value;
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
