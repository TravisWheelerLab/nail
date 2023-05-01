use crate::align::bounded::structs::{CloudBoundGroup, RowBounds};
use serde_json::json;

pub type JsonVec = Vec<serde_json::Value>;

impl CloudBoundGroup {
    pub fn json(&self, id: &str) -> (JsonVec, JsonVec) {
        let mut bounds: JsonVec = vec![];
        let mut diagonal: JsonVec = vec![];

        for bound in self.bounds() {
            bounds.push(json!({
                "id": format!("left-{}-{}", id, bound.anti_diagonal_idx()),
                "start": bound.left_profile_idx,
                "end": bound.left_profile_idx + 1,
                "row": bound.left_target_idx,
            }));
            bounds.push(json!({
                "id": format!("right-{}-{}", id, bound.anti_diagonal_idx()),
                "start": bound.right_profile_idx,
                "end": bound.right_profile_idx + 1,
                "row": bound.right_target_idx,
            }));
            diagonal.push(json!({
                    "id": format!("anti-diagonal-{}-{}", id, bound.anti_diagonal_idx()),
                    "start": 0,
                    "end": 0,
                    "x1" : bound.left_profile_idx as f32 + 0.5,
                    "x2" : bound.right_profile_idx as f32 + 0.5,
                    "y1" : bound.left_target_idx as f32 + 0.5,
                    "y2" : bound.right_target_idx as f32 + 0.5,
            }));
        }
        (bounds, diagonal)
    }
}

impl RowBounds {
    pub fn json(&self) -> JsonVec {
        let mut vec: JsonVec = vec![];
        for row in self.target_start..=self.target_end {
            vec.push(json!({
                "id": format!("row-{}", row),
                "start": self.left_row_bounds[row],
                "end": self.right_row_bounds[row] + 1,
                "row": row,
            }));
        }
        vec
    }
}
