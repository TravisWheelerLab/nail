use std::io::Write;

use crate::structs::profile::constants::{NUM_SPECIAL_STATES, SPECIAL_STATE_IDX_TO_NAME};
use anyhow::Result;

#[derive(Default)]
pub struct DpMatrix {
    pub profile_length: usize,
    pub target_length: usize,
    pub insert_matrix: Vec<Vec<f32>>,
    pub match_matrix: Vec<Vec<f32>>,
    pub delete_matrix: Vec<Vec<f32>>,
    pub special_matrix: Vec<Vec<f32>>,
}

impl DpMatrix {
    pub fn new(profile_length: usize, target_length: usize) -> Self {
        let mut matrix = DpMatrix {
            profile_length,
            target_length,
            insert_matrix: vec![vec![0.0; target_length + 1]; profile_length + 1],
            match_matrix: vec![vec![0.0; target_length + 1]; profile_length + 1],
            delete_matrix: vec![vec![0.0; target_length + 1]; profile_length + 1],
            special_matrix: vec![vec![0.0; 5]; profile_length + 1],
        };

        for i in 0..=target_length {
            matrix.set_match(i, 0, -f32::INFINITY);
            matrix.set_insert(i, 0, -f32::INFINITY);
            matrix.set_insert(i, profile_length, -f32::INFINITY);
            matrix.set_delete(i, 0, -f32::INFINITY);
            matrix.set_delete(i, 1, -f32::INFINITY);
        }

        matrix
    }

    #[inline(always)]
    pub fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.insert_matrix[target_idx][profile_idx]
    }

    #[inline(always)]
    pub fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.insert_matrix[target_idx][profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.match_matrix[target_idx][profile_idx]
    }

    #[inline(always)]
    pub fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.match_matrix[target_idx][profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.delete_matrix[target_idx][profile_idx]
    }

    #[inline(always)]
    pub fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.delete_matrix[target_idx][profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        self.special_matrix[target_idx][special_idx]
    }

    #[inline(always)]
    pub fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        self.special_matrix[target_idx][special_idx] = value;
    }

    pub fn dump(&self, out: &mut impl Write) -> Result<()> {
        let target_idx_width = self.target_length.to_string().len();
        let first_column_width = target_idx_width + 3;
        // TODO: should these be global statics or something?
        let column_width = 13;
        let precision = 8;

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
            write!(out, "{:.w$} M ", target_idx, w = target_idx_width)?;
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
