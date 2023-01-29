use std::cmp::min;

#[derive(Default)]
pub struct CloudAntiDiagonal {
    pub match_vector: Vec<f32>,
    pub insert_vector: Vec<f32>,
    pub delete_vector: Vec<f32>,
}

impl CloudAntiDiagonal {
    pub fn new(length: usize) -> Self {
        CloudAntiDiagonal {
            match_vector: vec![-f32::INFINITY; length + 2],
            insert_vector: vec![-f32::INFINITY; length + 2],
            delete_vector: vec![-f32::INFINITY; length + 2],
        }
    }
}

#[derive(Default)]
pub struct CloudMatrix {
    pub data: [CloudAntiDiagonal; 3],
}

impl CloudMatrix {
    pub fn new(profile_length: usize, target_length: usize) -> Self {
        let max_anti_diagonal_length = min(profile_length, target_length);
        CloudMatrix {
            data: [
                CloudAntiDiagonal::new(max_anti_diagonal_length),
                CloudAntiDiagonal::new(max_anti_diagonal_length),
                CloudAntiDiagonal::new(max_anti_diagonal_length),
            ],
        }
    }

    #[inline(always)]
    pub fn set_match(&mut self, row_idx: usize, col_idx: usize, value: f32) {
        self.data[row_idx].match_vector[col_idx] = value;
    }

    #[inline(always)]
    pub fn get_match(&self, row_idx: usize, col_idx: usize) -> f32 {
        self.data[row_idx].match_vector[col_idx]
    }

    #[inline(always)]
    pub fn set_insert(&mut self, row_idx: usize, col_idx: usize, value: f32) {
        self.data[row_idx].insert_vector[col_idx] = value;
    }

    #[inline(always)]
    pub fn get_insert(&self, row_idx: usize, col_idx: usize) -> f32 {
        self.data[row_idx].insert_vector[col_idx]
    }

    #[inline(always)]
    pub fn set_delete(&mut self, row_idx: usize, col_idx: usize, value: f32) {
        self.data[row_idx].delete_vector[col_idx] = value;
    }

    #[inline(always)]
    pub fn get_delete(&self, row_idx: usize, col_idx: usize) -> f32 {
        self.data[row_idx].delete_vector[col_idx]
    }
}
