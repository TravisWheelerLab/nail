use crate::align::bounded::structs::cloud_matrix::CloudAntiDiagonal;

#[derive(Default)]
pub struct CloudMatrixQuadratic {
    pub profile_length: usize,
    pub target_length: usize,
    pub row_data: Vec<CloudAntiDiagonal>,
    pub special_state_data: Vec<Vec<f32>>,
}
