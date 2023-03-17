#[derive(Default, Clone)]
pub struct CloudAntiDiagonal {
    pub match_vector: Vec<f32>,
    pub insert_vector: Vec<f32>,
    pub delete_vector: Vec<f32>,
}

impl CloudAntiDiagonal {
    pub fn new(length: usize) -> Self {
        CloudAntiDiagonal {
            // we buffer the length with +2 to pad each side
            match_vector: vec![-f32::INFINITY; length],
            insert_vector: vec![-f32::INFINITY; length],
            delete_vector: vec![-f32::INFINITY; length],
        }
    }

    pub fn resize(&mut self, new_length: usize) {
        self.match_vector.resize(new_length, -f32::INFINITY);
        self.insert_vector.resize(new_length, -f32::INFINITY);
        self.delete_vector.resize(new_length, -f32::INFINITY);
    }

    pub fn reset(&mut self) {
        for profile_idx in 0..self.match_vector.len() {
            self.match_vector[profile_idx] = -f32::INFINITY;
            self.insert_vector[profile_idx] = -f32::INFINITY;
            self.delete_vector[profile_idx] = -f32::INFINITY;
        }
    }
}

#[derive(Default)]
pub struct CloudMatrixLinear {
    profile_length: usize,
    col_count: usize,
    pub data: [CloudAntiDiagonal; 3],
}

impl CloudMatrixLinear {
    pub fn new(profile_length: usize) -> Self {
        let col_count = profile_length + 2;
        CloudMatrixLinear {
            profile_length,
            col_count,
            data: [
                CloudAntiDiagonal::new(col_count),
                CloudAntiDiagonal::new(col_count),
                CloudAntiDiagonal::new(col_count),
            ],
        }
    }

    pub fn resize(&mut self, new_profile_length: usize) {
        self.data[0].resize(new_profile_length + 2);
        self.data[1].resize(new_profile_length + 2);
        self.data[2].resize(new_profile_length + 2);
    }

    pub fn reset(&mut self) {
        self.data[0].reset();
        self.data[1].reset();
        self.data[2].reset();
    }

    pub fn reuse(&mut self, new_profile_length: usize) {
        let new_col_count = new_profile_length + 2;
        if new_col_count > self.col_count {
            self.resize(new_profile_length);
            self.col_count = new_col_count;
        }

        self.profile_length = new_profile_length;

        self.reset();
    }

    pub fn print(&self) {
        for (i, row) in self.data.iter().enumerate() {
            print!("M {}", i);
            for val in row.match_vector.iter() {
                print!("{:8.3} ", val)
            }
            println!();
            print!("I {}", i);
            for val in row.insert_vector.iter() {
                print!("{:8.3} ", val)
            }
            println!();
            print!("D {}", i);
            for val in row.delete_vector.iter() {
                print!("{:8.3} ", val)
            }
            println!();
            println!();
        }
        println!("----");
    }

    #[inline]
    pub fn set_match(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].match_vector[profile_idx] = value;
    }

    #[inline]
    pub fn get_match(&self, row_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].match_vector[profile_idx]
    }

    #[inline]
    pub fn set_insert(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].insert_vector[profile_idx] = value;
    }

    #[inline]
    pub fn get_insert(&self, row_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].insert_vector[profile_idx]
    }

    #[inline]
    pub fn set_delete(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].delete_vector[profile_idx] = value;
    }

    #[inline]
    pub fn get_delete(&self, row_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(profile_idx <= self.profile_length);
        self.data[row_idx].delete_vector[profile_idx]
    }
}
