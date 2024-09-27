use std::io::Write;

use crate::structs::Profile;
use anyhow::Result;

use super::RowBounds;

pub trait DpMatrix {
    fn target_length(&self) -> usize;
    fn profile_length(&self) -> usize;
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_special(&self, target_idx: usize, special_idx: usize) -> f32;
    fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32);
    fn dump(&self, out: &mut impl Write) -> Result<()> {
        let target_idx_width = self.target_length().to_string().len();
        let first_column_width = target_idx_width + 3;
        // TODO: should these be global statics or something?
        let column_width = 13;
        let precision = 3;

        // write the profile indices
        write!(out, "{}", " ".repeat(first_column_width - 1))?;
        for profile_idx in 0..=self.profile_length() {
            write!(out, "{:w$} ", profile_idx, w = column_width)?;
        }

        for special_idx in 0..Profile::NUM_SPECIAL_STATES {
            write!(
                out,
                "{:.w$} ",
                Profile::SPECIAL_STATE_IDX_TO_NAME[special_idx],
                w = column_width
            )?;
        }
        writeln!(out)?;

        write!(out, "{}", " ".repeat(first_column_width))?;
        for _ in 0..=self.profile_length() + Profile::NUM_SPECIAL_STATES {
            write!(out, "   {} ", "-".repeat(column_width - 3))?;
        }
        writeln!(out)?;

        for target_idx in 0..=self.target_length() {
            // write the match line
            write!(out, "{:w$} M ", target_idx, w = target_idx_width)?;
            for profile_idx in 0..=self.profile_length() {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_match(target_idx, profile_idx),
                    w = column_width,
                    p = precision
                )?;
            }

            // write the special states on the match line
            for special_idx in 0..Profile::NUM_SPECIAL_STATES {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_special(target_idx, special_idx),
                    w = column_width,
                    p = precision
                )?;
            }
            writeln!(out)?;

            // write the insert line
            write!(out, "{:w$} I ", target_idx, w = target_idx_width)?;
            for profile_idx in 0..=self.profile_length() {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_insert(target_idx, profile_idx),
                    w = column_width,
                    p = precision
                )?;
            }
            writeln!(out)?;

            // write the delete line
            write!(out, "{:w$} D ", target_idx, w = target_idx_width)?;
            for profile_idx in 0..=self.profile_length() {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_delete(target_idx, profile_idx),
                    w = column_width,
                    p = precision
                )?;
            }
            writeln!(out, "\n")?;
        }

        Ok(())
    }
}

#[derive(Default, Clone)]
pub struct DpMatrixFlat {
    pub target_length: usize,
    pub profile_length: usize,
    /// The DP matrix core model data cells as a flat vector.
    //
    // the data is stored in the following pattern:
    //     [
    //
    //         m_(0, 0), i_(0, 0), d_(0, 0),
    //         m_(0, 1), i_(0, 1), d_(0, 1),
    //         ...
    //         m_(0, P), i_(0, P), d_(0, P),
    //         ...
    //         m_(T, 0), i_(T, 0), d_(T, 0),
    //         ...
    //         m_(T, P), i_(T, P), d_(T, P)
    //
    //     ]
    //
    // where:
    //
    //     T:        <target_length>
    //     S:        <profile_length>
    //     m_(i, j): the match score at cell (i, j)
    //     i_(i, j): the insert score at cell (i, j)
    //     d_(i, j): the delete score at cell (i, j)
    //
    pub core_data: Vec<f32>,
    /// The DP matrix special state data cells as a flat vector.
    //
    // the data is stored in the following pattern:
    //     [
    //
    //         N_0, B_0, E_0, C_0, J_0,
    //         N_1, B_1, E_1, C_1, J_1,
    //         ...
    //         N_T, B_T, E_T, C_T, J_T,
    //
    //     ]
    //
    // where:
    //
    //     T:        <target_length>
    //     N_i: the N state score at target position i
    //     B_i: the B state score at target position i
    //     E_i: the E state score at target position i
    //     C_i: the C state score at target position i
    //     J_i: the J state score at target position i
    //
    pub special_data: Vec<f32>,
}

impl DpMatrixFlat {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        let core_length = 3 * (target_length + 1) * (profile_length + 1);
        let special_length = 5 * (target_length + 1);
        DpMatrixFlat {
            target_length,
            profile_length,
            core_data: vec![-f32::INFINITY; core_length],
            special_data: vec![-f32::INFINITY; special_length],
        }
    }

    pub fn resize(&mut self, new_target_length: usize, new_profile_length: usize) {
        let new_size = new_target_length * new_profile_length;
        if new_size > self.core_data.len() {
            self.core_data.resize(new_size, -f32::INFINITY);
        }
        self.target_length = new_target_length;
        self.profile_length = new_profile_length;
    }

    pub fn reset(&mut self) {
        for core_idx in 0..(3 * (self.target_length + 1) * (self.profile_length + 1)) {
            self.core_data[core_idx] = -f32::INFINITY;
        }

        for special_idx in 0..((self.target_length + 1) * 5) {
            self.special_data[special_idx] = -f32::INFINITY;
        }
    }

    pub fn reuse(&mut self, new_target_length: usize, new_profile_length: usize) {
        // TODO: need logic for resizing self.special_data
        let new_core_length = 3 * (new_target_length + 1) * (new_profile_length + 1);
        let new_special_length = 5 * (new_target_length + 1);

        if new_core_length > self.core_data.len() {
            panic!("tried to resize DpMatrix: core");
        } else if new_special_length > self.special_data.len() {
            panic!("tried to resize DpMatrix: special");
        }

        self.target_length = new_target_length;
        self.profile_length = new_profile_length;

        self.reset();
    }
}

impl DpMatrix for DpMatrixFlat {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (3 * profile_idx)]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (3 * profile_idx)] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 1)]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 1)] =
            value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 2)]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 2)] =
            value;
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

        self.core_data.resize(core_length, -f32::INFINITY);

        self.special_data.resize(new_special_length, -f32::INFINITY);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dp_matrix_flat() {
        let mut m = [[-f32::INFINITY; 6]; 6];
        let mut i = [[-f32::INFINITY; 6]; 6];
        let mut d = [[-f32::INFINITY; 6]; 6];

        (1..=5).for_each(|row| {
            (1..=5).for_each(|col| {
                m[row][col] = (row * 10 + col) as f32 / 10.0;
                i[row][col] = ((row + 5) * 10 + col) as f32 / 10.0;
                d[row][col] = ((row + 10) * 10 + col) as f32 / 10.0;
            });
        });

        let mut matrix = DpMatrixFlat::new(5, 5);

        (1..=5).for_each(|row| {
            (1..=5).for_each(|col| {
                matrix.set_match(row, col, m[row][col]);
                matrix.set_insert(row, col, i[row][col]);
                matrix.set_delete(row, col, d[row][col]);
            });
        });

        (0..=5).for_each(|row| {
            (0..=5).for_each(|col| {
                assert_eq!(matrix.get_match(row, col), m[row][col]);
                assert_eq!(matrix.get_insert(row, col), i[row][col]);
                assert_eq!(matrix.get_delete(row, col), d[row][col]);
            });
        });
    }

    #[test]
    fn test_dp_matrix_sparse() {
        let mut m = [[-f32::INFINITY; 6]; 6];
        let mut i = [[-f32::INFINITY; 6]; 6];
        let mut d = [[-f32::INFINITY; 6]; 6];

        let bounds = RowBounds {
            target_length: 1,
            profile_length: 5,
            target_start: 1,
            target_end: 5,
            row_capacity: 0,
            left_row_bounds: vec![0, 1, 1, 2, 3, 4],
            right_row_bounds: vec![0, 2, 3, 4, 5, 5],
            num_cells: 0,
        };

        (bounds.target_start..=bounds.target_end).for_each(|row| {
            (bounds.left_row_bounds[row]..=bounds.right_row_bounds[row]).for_each(|col| {
                m[row][col] = (row * 10 + col) as f32 / 10.0;
                i[row][col] = ((row + 5) * 10 + col) as f32 / 10.0;
                d[row][col] = ((row + 10) * 10 + col) as f32 / 10.0;
            });
        });

        let mut matrix = DpMatrixSparse::new(5, 5, &bounds);

        (bounds.target_start..=bounds.target_end).for_each(|row| {
            (bounds.left_row_bounds[row]..=bounds.right_row_bounds[row]).for_each(|col| {
                matrix.set_match(row, col, m[row][col]);
                matrix.set_insert(row, col, i[row][col]);
                matrix.set_delete(row, col, d[row][col]);
            });
        });

        (0..=5).for_each(|row| {
            (0..=5).for_each(|col| {
                assert_eq!(matrix.get_match(row, col), m[row][col]);
                assert_eq!(matrix.get_insert(row, col), i[row][col]);
                assert_eq!(matrix.get_delete(row, col), d[row][col]);
            });
        });
    }
}
