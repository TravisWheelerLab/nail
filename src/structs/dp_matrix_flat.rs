use std::io::Write;

use crate::structs::profile::constants::{NUM_SPECIAL_STATES, SPECIAL_STATE_IDX_TO_NAME};
use anyhow::Result;

#[derive(Default)]
pub struct DpMatrixFlat {
    pub target_length: usize,
    pub profile_length: usize,
    /// The DP matrix core model data cells as a flat vector.
    ///
    /// It's stored in the following pattern:
    ///     [
    ///
    ///         m_(0, 0), i_(0, 0), d_(0, 0),
    ///         m_(0, 1), i_(0, 1), d_(0, 1),
    ///         ...
    ///         m_(0, P), i_(0, P), d_(0, P),
    ///         ...
    ///         m_(T, 0), i_(T, 0), d_(T, 0),
    ///         ...
    ///         m_(T, P), i_(T, P), d_(T, P)
    ///
    ///     ]
    ///
    /// where:
    ///
    ///     T:        <target_length>
    ///     S:        <profile_length>
    ///     m_(i, j): the match score at cell (i, j)
    ///     i_(i, j): the insert score at cell (i, j)
    ///     d_(i, j): the delete score at cell (i, j)
    ///
    pub core_data: Vec<f32>,
    /// The DP matrix special state data cells as a flat vector.
    ///
    /// It's stored in the following pattern:
    ///     [
    ///
    ///         N_0, B_0, E_0, C_0, J_0,
    ///         N_1, B_1, E_1, C_1, J_1,
    ///         ...
    ///         N_T, B_T, E_T, C_T, J_T,
    ///
    ///     ]
    ///
    /// where:
    ///
    ///     T:        <target_length>
    ///     N_i: the N state score at target position i
    ///     B_i: the B state score at target position i
    ///     E_i: the E state score at target position i
    ///     C_i: the C state score at target position i
    ///     J_i: the J state score at target position i
    ///
    pub special_data: Vec<f32>,
}

impl DpMatrixFlat {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        DpMatrixFlat {
            target_length,
            profile_length,
            core_data: vec![],
            special_data: vec![],
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

    #[inline]
    pub fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + profile_idx]
    }

    #[inline]
    pub fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + profile_idx] = value;
    }

    #[inline]
    pub fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (profile_idx + 1)]
    }

    #[inline]
    pub fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (profile_idx + 1)] = value;
    }

    #[inline]
    pub fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (profile_idx + 2)]
    }

    #[inline]
    pub fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (profile_idx + 2)] = value;
    }

    #[inline]
    pub fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx]
    }

    #[inline]
    pub fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx] = value;
    }

    pub fn dump(&self, out: &mut impl Write) -> Result<()> {
        let target_idx_width = self.target_length.to_string().len();
        let first_column_width = target_idx_width + 3;
        // TODO: should these be global statics or something?
        let column_width = 13;
        let precision = 3;

        // write the profile indices
        write!(out, "{}", " ".repeat(first_column_width - 1))?;
        for profile_idx in 0..=self.profile_length {
            write!(out, "{:w$} ", profile_idx, w = column_width)?;
        }

        for special_idx in 0..NUM_SPECIAL_STATES {
            write!(
                out,
                "{:.w$} ",
                SPECIAL_STATE_IDX_TO_NAME[special_idx],
                w = column_width
            )?;
        }
        writeln!(out)?;

        write!(out, "{}", " ".repeat(first_column_width))?;
        for _ in 0..=self.profile_length + NUM_SPECIAL_STATES {
            write!(out, "   {} ", "-".repeat(column_width - 3))?;
        }
        writeln!(out)?;

        for target_idx in 0..=self.target_length {
            // write the match line
            write!(out, "{:w$} M ", target_idx, w = target_idx_width)?;
            for profile_idx in 0..=self.profile_length {
                write!(
                    out,
                    "{:w$.p$} ",
                    self.get_match(target_idx, profile_idx),
                    w = column_width,
                    p = precision
                )?;
            }

            // write the special states on the match line
            for special_idx in 0..NUM_SPECIAL_STATES {
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
            for profile_idx in 0..=self.profile_length {
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
            for profile_idx in 0..=self.profile_length {
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
