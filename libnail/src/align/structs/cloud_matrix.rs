use std::io::Write;

pub trait CloudMatrix {
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32;
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32);
    fn reset_anti_diagonal(&mut self, anti_diagonal_idx: usize);
    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()>;
}

pub struct CloudMatrixQuadratic {
    pub target_length: usize,
    pub profile_length: usize,
    pub match_data: Vec<Vec<f32>>,
    pub insert_data: Vec<Vec<f32>>,
    pub delete_data: Vec<Vec<f32>>,
}

impl CloudMatrixQuadratic {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        Self {
            target_length,
            profile_length,
            match_data: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
            insert_data: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
            delete_data: vec![vec![-f32::INFINITY; profile_length + 1]; target_length + 1],
        }
    }
}

impl CloudMatrixQuadratic {
    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
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

impl CloudMatrix for CloudMatrixQuadratic {
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.match_data[target_idx][profile_idx]
    }

    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.match_data[target_idx][profile_idx] = value;
    }

    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.insert_data[target_idx][profile_idx]
    }

    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.insert_data[target_idx][profile_idx] = value;
    }

    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.delete_data[target_idx][profile_idx]
    }

    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.delete_data[target_idx][profile_idx] = value;
    }

    fn reset_anti_diagonal(&mut self, _anti_diagonal_idx: usize) {
        // no-op for the quadratic matrix
    }

    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
        self.dump(out)
    }
}

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

#[derive(Default, Clone)]
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

    pub fn row_idx(target_idx: usize, profile_idx: usize) -> usize {
        (target_idx + profile_idx) % 3
    }
}

impl CloudMatrix for CloudMatrixLinear {
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.get_match(Self::row_idx(target_idx, profile_idx), profile_idx)
    }

    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.set_match(Self::row_idx(target_idx, profile_idx), profile_idx, value);
    }

    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.get_insert(Self::row_idx(target_idx, profile_idx), profile_idx)
    }

    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.set_insert(Self::row_idx(target_idx, profile_idx), profile_idx, value);
    }

    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        self.get_delete(Self::row_idx(target_idx, profile_idx), profile_idx)
    }

    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        self.set_delete(Self::row_idx(target_idx, profile_idx), profile_idx, value);
    }

    fn reset_anti_diagonal(&mut self, anti_diagonal_idx: usize) {
        self.data[anti_diagonal_idx].reset();
    }

    fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
        panic!("can't dump a linear cloud matrix");
    }
}
