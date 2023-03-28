use std::io::Write;

use crate::structs::Profile;
use anyhow::Result;

pub trait DpMatrix {
    fn target_length(&self) -> usize;
    fn profile_length(&self) -> usize;
    fn resize(&mut self, new_target_length: usize, new_profile_length: usize);
    fn reset(&mut self);
    fn reuse(&mut self, new_target_length: usize, new_profile_length: usize);
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_special(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_special(&mut self, target_idx: usize, profile_idx: usize, value: f32);
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

#[derive(Default)]
pub struct DpMatrix3D {
    pub target_length: usize,
    pub profile_length: usize,
    pub row_count: usize,
    pub col_count: usize,
    pub insert_matrix: Vec<Vec<f32>>,
    pub match_matrix: Vec<Vec<f32>>,
    pub delete_matrix: Vec<Vec<f32>>,
    pub special_matrix: Vec<Vec<f32>>,
}

impl DpMatrix3D {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        DpMatrix3D {
            target_length,
            profile_length,
            row_count: target_length + 1,
            col_count: profile_length + 1,
            insert_matrix: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
            match_matrix: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
            delete_matrix: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
            special_matrix: vec![vec![-f32::INFINITY; 5]; target_length + 1],
        }
    }
}

impl DpMatrix for DpMatrix3D {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    fn resize(&mut self, new_target_length: usize, new_profile_length: usize) {
        #![allow(unused_variables)]
        todo!()
    }

    fn reset(&mut self) {
        for target_idx in 0..=self.target_length {
            for special_state_idx in 0..Profile::NUM_SPECIAL_STATES {
                self.set_special(target_idx, special_state_idx, -f32::INFINITY);
            }
            for profile_idx in 0..=self.profile_length {
                self.set_match(target_idx, profile_idx, -f32::INFINITY);
                self.set_insert(target_idx, profile_idx, -f32::INFINITY);
                self.set_delete(target_idx, profile_idx, -f32::INFINITY);
            }
        }
    }

    fn reuse(&mut self, new_target_length: usize, new_profile_length: usize) {
        let new_row_count = new_target_length + 1;
        let new_col_count = new_profile_length + 1;

        if new_row_count > self.row_count || new_col_count > self.col_count {
            panic!("tried to resize DpMatrix");
        }

        self.target_length = new_target_length;
        self.profile_length = new_profile_length;

        self.reset();
    }

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.match_matrix[target_idx][profile_idx]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.match_matrix[target_idx][profile_idx] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.insert_matrix[target_idx][profile_idx]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.insert_matrix[target_idx][profile_idx] = value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.delete_matrix[target_idx][profile_idx]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.delete_matrix[target_idx][profile_idx] = value;
    }

    #[inline]
    fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_matrix[target_idx][special_idx]
    }

    #[inline]
    fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_matrix[target_idx][special_idx] = value;
    }
}
