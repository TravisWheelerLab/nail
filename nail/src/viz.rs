use std::fs;
use std::path::Path;

use libnail::align::structs::{CloudBoundGroup, RowBounds};
use serde::Serialize;

#[allow(dead_code)]
pub fn write_soda_html(
    data: &impl Serialize,
    template_path: impl AsRef<Path>,
    js_path: impl AsRef<Path>,
    out_path: impl AsRef<Path>,
) {
    let template = fs::read_to_string(template_path).expect("failed to read template");
    let js = fs::read_to_string(js_path).expect("failed to read js");

    let mut viz_html = template.replace(
        "DATA_TARGET",
        &serde_json::to_string(data).expect("failed to serialize JSON data"),
    );

    viz_html = viz_html.replace("JS_TARGET", &js);

    let mut file = std::fs::File::create(out_path).expect("failed to create file");

    std::io::Write::write_all(&mut file, viz_html.as_bytes()).expect("failed to write to file");
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AntiDiagonalBoundSodaData {
    row_start: usize,
    row_end: usize,
    col_start: usize,
    col_end: usize,
    forward_bounds: (String, String),
    backward_bounds: (String, String),
}

#[allow(dead_code)]
impl AntiDiagonalBoundSodaData {
    pub fn new(forward: &CloudBoundGroup, backward: &CloudBoundGroup) -> Self {
        let row_start = forward
            .get_first()
            .right_target_idx
            .min(backward.get_first().right_target_idx);

        let col_start = forward
            .get_first()
            .left_profile_idx
            .min(backward.get_first().left_profile_idx);

        let row_end = forward
            .get_last()
            .left_target_idx
            .max(backward.get_last().left_target_idx);

        let col_end = forward
            .get_last()
            .right_profile_idx
            .max(backward.get_last().right_profile_idx);

        let mut forward_bounds: (String, String) = forward
            .bounds()
            .iter()
            .map(|b| {
                (
                    format!("{},{}|", b.left_target_idx, b.left_profile_idx,),
                    format!("{},{}|", b.right_target_idx, b.right_profile_idx),
                )
            })
            .unzip();

        let mut backward_bounds: (String, String) = backward
            .bounds()
            .iter()
            .map(|b| {
                (
                    format!("{},{}|", b.left_target_idx, b.left_profile_idx,),
                    format!("{},{}|", b.right_target_idx, b.right_profile_idx),
                )
            })
            .unzip();

        forward_bounds.0.pop();
        forward_bounds.1.pop();
        backward_bounds.0.pop();
        backward_bounds.1.pop();

        Self {
            row_start,
            row_end,
            col_start,
            col_end,
            forward_bounds,
            backward_bounds,
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RowBoundSodaData {
    row_start: usize,
    row_end: usize,
    col_start: usize,
    col_end: usize,
    row_bounds: String,
}

#[allow(dead_code)]
impl RowBoundSodaData {
    pub fn new(bounds: &RowBounds) -> Self {
        Self {
            row_start: todo!(),
            row_end: todo!(),
            col_start: todo!(),
            col_end: todo!(),
            row_bounds: todo!(),
        }
    }
}
