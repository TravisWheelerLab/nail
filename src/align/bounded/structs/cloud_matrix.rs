use std::cmp::max;

#[derive(Default)]
pub struct CloudAntiDiagonal {
    pub match_vector: Vec<f32>,
    pub insert_vector: Vec<f32>,
    pub delete_vector: Vec<f32>,
}

impl CloudAntiDiagonal {
    pub fn new(length: usize) -> Self {
        CloudAntiDiagonal {
            // we buffer the length with +2 to pad each side
            match_vector: vec![-f32::INFINITY; length + 2],
            insert_vector: vec![-f32::INFINITY; length + 2],
            delete_vector: vec![-f32::INFINITY; length + 2],
        }
    }

    pub fn reuse(&mut self) {
        for i in 0..self.match_vector.len() {
            self.match_vector[i] = -f32::INFINITY;
            self.insert_vector[i] = -f32::INFINITY;
            self.delete_vector[i] = -f32::INFINITY;
        }
    }
}

#[derive(Default)]
pub struct CloudMatrix {
    pub data: [CloudAntiDiagonal; 3],
}

impl CloudMatrix {
    pub fn new(profile_length: usize, target_length: usize) -> Self {
        let length = max(profile_length, target_length);
        CloudMatrix {
            data: [
                CloudAntiDiagonal::new(length),
                CloudAntiDiagonal::new(length),
                CloudAntiDiagonal::new(length),
            ],
        }
    }

    pub fn reuse(&mut self) {
        self.data[0].reuse();
        self.data[1].reuse();
        self.data[2].reuse();
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

    #[inline(always)]
    pub fn set_match(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        self.data[row_idx].match_vector[profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_match(&self, row_idx: usize, profile_idx: usize) -> f32 {
        self.data[row_idx].match_vector[profile_idx]
    }

    #[inline(always)]
    pub fn set_insert(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        self.data[row_idx].insert_vector[profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_insert(&self, row_idx: usize, profile_idx: usize) -> f32 {
        self.data[row_idx].insert_vector[profile_idx]
    }

    #[inline(always)]
    pub fn set_delete(&mut self, row_idx: usize, profile_idx: usize, value: f32) {
        self.data[row_idx].delete_vector[profile_idx] = value;
    }

    #[inline(always)]
    pub fn get_delete(&self, row_idx: usize, profile_idx: usize) -> f32 {
        self.data[row_idx].delete_vector[profile_idx]
    }
}
