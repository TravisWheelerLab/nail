use std::{
    io::Write,
    ops::{Index, IndexMut},
};

use crate::{align::structs::AntiDiagonal, structs::Profile};

use super::{Cloud, RowBounds};

#[derive(Debug)]
#[repr(usize)]
pub enum CoreState {
    M(usize) = 0,
    I(usize) = 1,
    D(usize) = 2,
}

impl CoreState {
    #[inline(always)]
    fn state_idx(&self) -> usize {
        // safety:
        //   when an enum is repr(usize), it is guaranteed that
        //   the discriminant can be reliably cast from from
        //   the first 8 bytes of the in-memory representation
        //
        //   see RFC 2368: arbitrary_enum_discriminant:
        //     https://rust-lang.github.io/rfcs/2363-arbitrary-enum-discriminant.html
        //     https://doc.rust-lang.org/reference/items/enumerations.html#pointer-casting
        let raw_ptr = <*const _>::from(self).cast::<usize>();
        unsafe { *raw_ptr }
    }

    #[inline(always)]
    fn profile_idx(&self) -> usize {
        // safety:
        //   when an enum is repr(usize), it is guaranteed that
        //   the discriminant can be reliably cast from from
        //   the first 8 bytes of the in-memory representation
        //
        //   see RFC 2368: arbitrary_enum_discriminant:
        //     https://rust-lang.github.io/rfcs/2363-arbitrary-enum-discriminant.html
        //     https://doc.rust-lang.org/reference/items/enumerations.html#pointer-casting
        let mut raw_ptr = <*const _>::from(self).cast::<usize>();
        unsafe {
            raw_ptr = raw_ptr.add(1);
            *raw_ptr
        }
    }
}

#[repr(usize)]
pub enum BackgroundState {
    N = 0,
    B = 1,
    C = 2,
    E = 3,
}

pub type CoreCell = (CoreState, usize);
pub type BackgroundCell = (BackgroundState, usize);

pub trait CoreCellIndexable:
    Index<CoreCell, Output = f32> + IndexMut<CoreCell, Output = f32>
{
}
pub trait BackgroundCellIndexable:
    Index<BackgroundCell, Output = f32> + IndexMut<BackgroundCell, Output = f32>
{
}
pub trait CellIndexable: CoreCellIndexable + BackgroundCellIndexable {}

#[derive(Default, Clone)]
pub struct AdMatrixSparse {
    target_len: usize,
    profile_len: usize,
    first_ad_idx: usize,
    block_offsets: Vec<usize>,
    ad_start_offsets: Vec<usize>,
    core_data: [Vec<f32>; 3],
    background_data: [Vec<f32>; 4],
}

impl AdMatrixSparse {
    pub fn layout(&self) {
        println!("target len:     {}", self.target_len());
        println!("profile len:    {}", self.profile_len());
        println!("first AD index: {}", self.first_ad_idx);

        self.block_offsets
            .iter()
            .zip(self.block_offsets.iter().skip(1))
            .zip(self.ad_start_offsets.iter())
            .enumerate()
            .for_each(|(idx, ((&a, &b), o))| {
                print!("AD: {:<3} -> {a:<3} | +{o:<3} | ", idx + self.first_ad_idx);
                self.core_data[0][(a + 1)..b].iter().for_each(|v| {
                    print!("{v:4.2} ");
                });
                println!();
            });
        println!(
            "AD: ∅   -> {:<3} | +{:<3} | {:4.2}",
            self.block_offsets.last().unwrap(),
            0,
            -f32::INFINITY
        );
    }

    pub fn from_cloud(cloud: &Cloud) -> Self {
        let mut matrix = Self::default();
        matrix.reuse(cloud);
        matrix
    }

    pub fn reuse(&mut self, cloud: &Cloud) {
        self.target_len = cloud.target_len();
        self.profile_len = cloud.profile_len();

        self.first_ad_idx = cloud.first().idx();

        let num_ad = cloud.num_anti_diagonals();
        // +1 to point to the final pad cell
        let num_blocks = num_ad + 1;
        self.ad_start_offsets.resize(num_ad, 0);
        self.block_offsets.resize(num_blocks, 0);
        self.ad_start_offsets.shrink_to_fit();
        self.block_offsets.shrink_to_fit();

        /// A simple local function to compute the AD start offset
        /// given a `Bound` and the length of a `Profile`
        #[inline(always)]
        fn offset(bound: &AntiDiagonal, profile_len: usize) -> usize {
            let left = bound.right_target_major();
            let delta = left.idx().saturating_sub(profile_len);
            println!("{bound:?}");
            println!("Left   {left:?}");
            println!("delta: {delta}\n");
            left.seq_idx - delta
        }

        // peeling off the first pair of offsets so
        // that the following iterator stays clean
        self.block_offsets[0] = 0;
        self.ad_start_offsets[0] = offset(cloud.first(), self.profile_len);

        cloud
            .bounds()
            .iter()
            .rev()
            .enumerate()
            .zip(cloud.bounds().iter().rev().enumerate().skip(1))
            .for_each(|((prev_idx, prev_bound), (idx, bound))| {
                self.ad_start_offsets[idx] = offset(bound, self.profile_len);
                let block_offset = self.block_offsets[prev_idx] + prev_bound.len() + 1;
                self.block_offsets[idx] = block_offset;
            });

        self.block_offsets[num_blocks - 1] =
            self.block_offsets[num_ad - 1] + cloud.last().len() + 1;

        let num_cells = self.block_offsets[num_blocks - 1] + 1;

        self.core_data.iter_mut().for_each(|v| {
            v.resize(num_cells, 0.0);
            v.shrink_to_fit();
        });

        self.background_data.iter_mut().for_each(|v| {
            v.resize(self.target_len + 1, 0.0);
            v.shrink_to_fit();
        });

        self.core_data
            .iter_mut()
            .for_each(|state_vec| state_vec.iter_mut().for_each(|val| *val = -f32::INFINITY));
    }

    pub fn size(&self) -> usize {
        self.core_data[0].len()
    }

    #[inline(always)]
    fn sparse_core_idx(&self, cell: &CoreCell) -> usize {
        let profile_idx = cell.0.profile_idx();
        let seq_idx = cell.1;

        let ad_idx = profile_idx + seq_idx - self.first_ad_idx;
        let seq_delta = ad_idx.saturating_sub(self.profile_len);
        let seq_idx_in_block = seq_idx - seq_delta;
        let block_start = self.block_offsets[ad_idx] + 1;
        let next_block_start = self.block_offsets[ad_idx + 1];
        let ad_offset = self.ad_start_offsets[ad_idx];

        // the saturating sub prevents moving left out of the block
        // the min with the next block offset prevents moving right out of the block
        (block_start + (seq_idx_in_block.saturating_sub(ad_offset))).min(next_block_start)
    }
}

impl Index<CoreCell> for AdMatrixSparse {
    type Output = f32;

    fn index(&self, cell: CoreCell) -> &Self::Output {
        let state_idx = cell.0.state_idx();
        let sparse_idx = self.sparse_core_idx(&cell);
        &self.core_data[state_idx][sparse_idx]
    }
}

impl IndexMut<CoreCell> for AdMatrixSparse {
    fn index_mut(&mut self, cell: CoreCell) -> &mut Self::Output {
        let state_idx = cell.0.state_idx();
        let sparse_idx = self.sparse_core_idx(&cell);
        &mut self.core_data[state_idx][sparse_idx]
    }
}

impl Index<BackgroundCell> for AdMatrixSparse {
    type Output = f32;

    fn index(&self, cell: BackgroundCell) -> &Self::Output {
        &self.background_data[cell.0 as usize][cell.1]
    }
}

impl IndexMut<BackgroundCell> for AdMatrixSparse {
    fn index_mut(&mut self, cell: BackgroundCell) -> &mut Self::Output {
        &mut self.background_data[cell.0 as usize][cell.1]
    }
}

impl CoreCellIndexable for AdMatrixSparse {}
impl BackgroundCellIndexable for AdMatrixSparse {}
impl CellIndexable for AdMatrixSparse {}

impl NewDpMatrix for AdMatrixSparse {
    fn target_len(&self) -> usize {
        self.target_len
    }

    fn profile_len(&self) -> usize {
        self.profile_len
    }
}

pub trait NewDpMatrix: CellIndexable {
    fn target_len(&self) -> usize;
    fn profile_len(&self) -> usize;

    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
        use BackgroundState::*;
        use CoreState::*;

        let t_idx_width = self.target_len().to_string().len();
        let first_column_width = t_idx_width + 3;
        // TODO: configurable?
        const W: usize = 13;
        const P: usize = 3;

        // -- profile indices
        write!(out, "{}", " ".repeat(first_column_width - 1))?;
        for p_idx in 0..=self.profile_len() {
            write!(out, "{:w$} ", p_idx, w = W)?;
        }

        for s_idx in 0..Profile::NUM_SPECIAL_STATES {
            write!(
                out,
                "{:>w$} ",
                Profile::SPECIAL_STATE_IDX_TO_NAME[s_idx],
                w = W
            )?;
        }
        writeln!(out)?;

        write!(out, "{}", " ".repeat(first_column_width))?;
        for _ in 0..=self.profile_len() + Profile::NUM_SPECIAL_STATES {
            write!(out, "   {} ", "-".repeat(W - 3))?;
        }
        writeln!(out)?;

        for t_idx in 0..=self.target_len() {
            // -- match
            write!(out, "{:w$} M ", t_idx, w = t_idx_width)?;
            for p_idx in 0..=self.profile_len() {
                write!(out, "{:w$.p$} ", self[(M(p_idx), t_idx)], w = W, p = P)?;
            }

            // -- special
            write!(out, "{:w$.p$} ", self[(E, t_idx)], w = W, p = P)?;
            write!(out, "{:w$.p$} ", self[(N, t_idx)], w = W, p = P)?;
            write!(out, "{:w$.p$} ", -f32::INFINITY, w = W, p = P)?;
            write!(out, "{:w$.p$} ", self[(B, t_idx)], w = W, p = P)?;
            write!(out, "{:w$.p$} ", self[(C, t_idx)], w = W, p = P)?;
            writeln!(out)?;

            // -- insert
            write!(out, "{:w$} I ", t_idx, w = t_idx_width)?;
            for p_idx in 0..=self.profile_len() {
                write!(out, "{:w$.p$} ", self[(I(p_idx), t_idx)], w = W, p = P)?;
            }
            writeln!(out)?;

            // -- delete
            write!(out, "{:w$} D ", t_idx, w = t_idx_width)?;
            for p_idx in 0..=self.profile_len() {
                write!(out, "{:w$.p$} ", self[(D(p_idx), t_idx)], w = W, p = P)?;
            }
            writeln!(out, "\n")?;
        }

        Ok(())
    }
}

impl<T: NewDpMatrix> DpMatrix for T {
    fn target_length(&self) -> usize {
        self.target_len()
    }

    fn profile_length(&self) -> usize {
        self.profile_len()
    }

    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self[(CoreState::M(profile_idx), target_idx)]
    }

    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self[(CoreState::M(profile_idx), target_idx)] = value;
    }

    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self[(CoreState::I(profile_idx), target_idx)]
    }

    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self[(CoreState::I(profile_idx), target_idx)] = value;
    }

    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self[(CoreState::D(profile_idx), target_idx)]
    }

    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self[(CoreState::D(profile_idx), target_idx)] = value;
    }

    fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        let state: BackgroundState = unsafe { std::mem::transmute(special_idx) };
        self[(state, target_idx)]
    }

    fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        let state: BackgroundState = unsafe { std::mem::transmute(special_idx) };
        self[(state, target_idx)] = value;
    }
}

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
    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
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
                "{:>w$} ",
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
            for special_idx in 0..(Profile::NUM_SPECIAL_STATES - 1) {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_special(target_idx, special_idx),
                    w = column_width,
                    p = precision
                )?;
            }
            write!(
                out,
                "{:w$.p$} ",
                -f32::INFINITY,
                w = column_width,
                p = precision
            )?;

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

        self.core_data.shrink_to_fit();
        self.special_data.shrink_to_fit();

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
    use super::CoreState::*;
    use super::*;
    use crate::align::structs::test_consts::BOUNDS_A_1;

    use assert2::assert;

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
                assert!(matrix.get_match(row, col) == m[row][col]);
                assert!(matrix.get_insert(row, col) == i[row][col]);
                assert!(matrix.get_delete(row, col) == d[row][col]);
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
                assert!(matrix.get_match(row, col) == m[row][col]);
                assert!(matrix.get_insert(row, col) == i[row][col]);
                assert!(matrix.get_delete(row, col) == d[row][col]);
            });
        });
    }

    #[test]
    fn test_ad_matrix_sparse_square() -> anyhow::Result<()> {
        const S: usize = 3;
        const P: usize = 3;

        let mut cloud = Cloud::new(S, P);
        cloud.fill_rectangle(0, 0, S, P);
        let mut mx = AdMatrixSparse::from_cloud(&cloud);

        let mut val = 1.0;
        (1..=P).for_each(|p| {
            (1..=S).for_each(|s| {
                mx[(M(p), s)] = val;
                mx[(I(p), s)] = val + 0.1;
                mx[(D(p), s)] = val + 0.2;
                val += 1.0;
            })
        });

        // the indices and respective values
        // that should be set by the loop above
        let set_idx = [7, 12, 17, 11, 16, 20, 19, 15, 22];
        let set_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 7.0, 9.0];

        set_idx.iter().zip(set_values).for_each(|(&idx, val)| {
            assert!(mx.core_data[0][idx] == val);
            assert!(mx.core_data[1][idx] == val + 0.1);
            assert!(mx.core_data[2][idx] == val + 0.2);
        });

        (0..mx.size())
            // everything that wasn't set should be -inf
            .filter(|idx| !set_idx.contains(idx))
            .for_each(|idx| {
                assert!(mx.core_data[0][idx] == -f32::INFINITY);
                assert!(mx.core_data[1][idx] == -f32::INFINITY);
                assert!(mx.core_data[2][idx] == -f32::INFINITY);
            });

        Ok(())
    }

    #[test]
    fn test_ad_matrix_sparse_cloud() -> anyhow::Result<()> {
        const S: usize = 10;
        const P: usize = 10;

        let mut cloud = Cloud::new(S, P);
        cloud.fill(&BOUNDS_A_1)?;
        cloud.image().print();
        let mut mx = AdMatrixSparse::from_cloud(&cloud);
        mx.layout();

        let mut val = 1.0;
        (1..=P).for_each(|p| {
            (1..=S).for_each(|s| {
                mx[(M(p), s)] = val;
                mx[(I(p), s)] = val + 0.1;
                mx[(D(p), s)] = val + 0.2;
                val += 1.0;
            })
        });

        let set_idx = [];
        let set_values = [];

        set_idx.iter().zip(set_values).for_each(|(&idx, val)| {
            assert!(mx.core_data[0][idx] == val);
            assert!(mx.core_data[1][idx] == val + 0.1);
            assert!(mx.core_data[2][idx] == val + 0.2);
        });

        (0..mx.size())
            .filter(|idx| !set_idx.contains(idx))
            .for_each(|idx| {
                assert!(mx.core_data[0][idx] == -f32::INFINITY);
                assert!(mx.core_data[1][idx] == -f32::INFINITY);
                assert!(mx.core_data[2][idx] == -f32::INFINITY);
            });

        Ok(())
    }
}
