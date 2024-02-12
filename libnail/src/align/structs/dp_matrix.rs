use std::io::Write;

use crate::structs::Profile;
use anyhow::Result;

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
