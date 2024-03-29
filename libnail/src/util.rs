use lazy_static::lazy_static;

pub trait Average<T> {
    fn avg(&self) -> T;
}

impl Average<usize> for Vec<usize> {
    fn avg(&self) -> usize {
        let sum: usize = self.iter().sum();
        sum / self.len()
    }
}

pub trait PrintMe {
    fn print(&self);
}

impl<T: PrintMe> PrintMe for Vec<T> {
    fn print(&self) {
        for val in self.iter() {
            val.print();
        }
    }
}

impl<T: PrintMe> PrintMe for &[T] {
    fn print(&self) {
        for val in self.iter() {
            val.print();
        }
    }
}

impl PrintMe for String {
    fn print(&self) {
        println!("{}", self)
    }
}

impl PrintMe for usize {
    fn print(&self) {
        println!("{}", self)
    }
}

impl PrintMe for i32 {
    fn print(&self) {
        println!("{}", self)
    }
}

impl PrintMe for f32 {
    fn print(&self) {
        println!("{:.3}", self)
    }
}

pub trait LogAbuse {
    fn ln_or_inf(self) -> f32;
    fn ln_or_max(self) -> f32;
}

impl LogAbuse for f32 {
    fn ln_or_inf(self) -> f32 {
        if self == 0.0 {
            -f32::INFINITY
        } else {
            self.ln()
        }
    }

    fn ln_or_max(self) -> f32 {
        if self == 0.0 {
            -f32::MAX
        } else {
            self.ln()
        }
    }
}

pub fn f32_vec_argmax(vec: &Vec<f32>) -> usize {
    let mut max: f32 = vec[0];
    let mut argmax: usize = 0;

    for i in 1..vec.len() {
        if vec[i] > max {
            max = vec[i];
            argmax = i;
        }
    }
    argmax
}

lazy_static! {
    pub static ref LOGSUM_LOOKUP: Vec<f32> = {
        let mut f: Vec<f32> = vec![];
        for i in 0..LOGSUM_TABLE_SIZE {
            f.push((1.0 + (-(i as f64) / LOGSUM_SCALE as f64).exp()).ln() as f32);
        }
        f
    };
}

const LOGSUM_SCALE: f32 = 1000.0;
const LOGSUM_TABLE_SIZE: usize = 16000;

/// A fast, table driven approximation of the sum of two floats in log space.
#[inline(always)]
pub fn log_add(a: f32, b: f32) -> f32 {
    let min = f32::min(a, b);
    let max = f32::max(a, b);

    debug_assert!(!a.is_nan());
    debug_assert!(!b.is_nan());
    debug_assert!(!a.is_sign_positive() || a.is_finite());
    debug_assert!(!b.is_sign_positive() || b.is_finite());

    if min == -f32::INFINITY || max - min >= 15.7 {
        max
    } else {
        max + LOGSUM_LOOKUP[((max - min) * LOGSUM_SCALE) as usize]
    }
}

#[macro_export]
macro_rules! log_sum {
    // Base case:
    ($x:expr) => ($x);
    // `$x` followed by at least one `$y,`
    ($x:expr, $($y:expr),+) => (
        // Call `log_sum!` on the tail `$y`
        log_add($x, log_sum!($($y),+))
    )
}

#[macro_export]
macro_rules! max_f32 {
    // Base case:
    ($x:expr) => ($x);
    // `$x` followed by at least one `$y,`
    ($x:expr, $($y:expr),+) => (
        // Call `max_f32!` on the tail `$y`
        $x.max(max_f32!($($y),+))
    )
}
