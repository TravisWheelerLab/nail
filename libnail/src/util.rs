use std::{
    cmp::Ordering::{Equal, Greater, Less},
    fmt::{Debug, Display},
    ops::{AddAssign, Div, DivAssign, MulAssign, SubAssign},
};

use lazy_static::lazy_static;

#[cfg(test)]
#[ctor::ctor]
fn init_backtrace() {
    color_backtrace::install();
}

pub trait Print {
    fn print(&self);
    fn print_debug(&self);
}

impl<T: Display + Debug> Print for T {
    fn print(&self) {
        println!("{self}");
    }

    fn print_debug(&self) {
        println!("{self:?}");
    }
}

pub trait CollectionPrint {
    fn print(&self);
    fn print_debug(&self);
}

impl<T: Display + Debug> CollectionPrint for Vec<T> {
    fn print(&self) {
        self.iter()
            .enumerate()
            .for_each(|(i, e)| println!("{i}: {e}"));
    }

    fn print_debug(&self) {
        self.iter()
            .enumerate()
            .for_each(|(i, e)| println!("{i}: {e:?}"));
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

pub trait Float:
    PartialOrd + Copy + Div<Output = Self> + AddAssign + SubAssign + DivAssign + MulAssign
{
    fn from_usize(n: usize) -> Self;
}

impl Float for f32 {
    fn from_usize(n: usize) -> Self {
        n as f32
    }
}

impl Float for f64 {
    fn from_usize(n: usize) -> Self {
        n as f64
    }
}

pub trait VecMath<T>
where
    T: Float,
{
    fn avg(&self) -> Option<T>;
    fn argmax(&self) -> Option<usize>;
    fn normalize(&mut self);
    fn add(&mut self, other: &[T]);
    fn sub(&mut self, other: &[T]);
    fn scale(&mut self, factor: T);
    fn saturate_lower(&mut self, min: T);
}

impl<T> VecMath<T> for Vec<T>
where
    T: Float,
{
    fn avg(&self) -> Option<T> {
        let mut sum = T::from_usize(0);
        self.iter().for_each(|&item| sum += item);

        Some(sum / T::from_usize(self.len()))
    }

    fn argmax(&self) -> Option<usize> {
        let mut max = *self.first()?;
        let mut argmax: usize = 0;

        for (idx, &item) in self.iter().enumerate().skip(1) {
            if item > max {
                max = item;
                argmax = idx;
            }
        }

        Some(argmax)
    }

    fn normalize(&mut self) {
        let mut sum = T::from_usize(0);
        self.iter().for_each(|&item| sum += item);
        self.iter_mut().for_each(|item| *item /= sum);
    }

    fn add(&mut self, other: &[T]) {
        self.iter_mut().zip(other).for_each(|(a, &b)| *a += b);
    }

    fn sub(&mut self, other: &[T]) {
        self.iter_mut().zip(other).for_each(|(a, &b)| *a -= b);
    }

    fn scale(&mut self, factor: T) {
        self.iter_mut().for_each(|item| *item *= factor);
    }

    fn saturate_lower(&mut self, min: T) {
        self.iter_mut().for_each(|item| {
            if *item < min {
                *item = min
            }
        });
    }
}

pub trait VecUtils<T>
where
    T: Clone,
{
    fn reset(&mut self, value: T);
    fn grow_or_shrink(&mut self, new_len: usize, value: T);
    fn resize_and_reset(&mut self, new_len: usize, value: T);
}

impl<T> VecUtils<T> for Vec<T>
where
    T: Clone,
{
    fn reset(&mut self, value: T) {
        self.iter_mut().for_each(|v| *v = value.clone());
    }

    fn resize_and_reset(&mut self, new_len: usize, value: T) {
        match new_len.cmp(&self.len()) {
            Less => {
                self.truncate(new_len);
                self.shrink_to_fit();
                self.iter_mut().for_each(|v| *v = value.clone());
            }
            Equal => self.iter_mut().for_each(|v| *v = value.clone()),
            Greater => {
                self.iter_mut().for_each(|v| *v = value.clone());
                self.resize(new_len, value);
            }
        }
    }

    fn grow_or_shrink(&mut self, new_len: usize, value: T) {
        match new_len.cmp(&self.len()) {
            Less => {
                self.truncate(new_len);
                self.shrink_to_fit();
            }
            Equal => {}
            Greater => {
                self.resize(new_len, value);
            }
        }
    }
}

pub fn mean_relative_entropy(a: &[Vec<f32>], b: &[f32]) -> f32 {
    a.iter()
        .map(|p| relative_entropy(p, b))
        .collect::<Vec<f32>>()
        .avg()
        .unwrap()
}

pub fn relative_entropy(a: &[f32], b: &[f32]) -> f32 {
    a.iter()
        .zip(b)
        .map(|(p_a, p_b)| p_a * (p_a / p_b).log2())
        .sum()
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

#[macro_export]
macro_rules! assert_eq_pairs {
    ($( $left:expr, $right:expr );+ $(;)?) => {
        $(
            assert!($left == $right);
        )+
    };
}

#[cfg(feature = "debug")]
pub mod debug {
    use std::{env, fs, path::PathBuf};

    use anyhow::anyhow;

    #[allow(dead_code)]
    pub fn workspace_root() -> anyhow::Result<PathBuf> {
        let mut dir = env::current_dir().unwrap();

        while dir.join("Cargo.toml").exists() {
            if dir.join("Cargo.toml").exists() {
                let content = fs::read_to_string(dir.join("Cargo.toml")).unwrap();
                if content.contains("[workspace]") {
                    return Ok(dir);
                }
            }
            dir = dir
                .parent()
                .ok_or(anyhow!("failed to find workspace root"))?
                .to_path_buf();
        }

        Err(anyhow!("failed to find workspace root"))
    }

    #[allow(dead_code)]
    pub fn debug_dir() -> anyhow::Result<PathBuf> {
        let dir = env::current_dir()?.join("debug/");

        if !dir.exists() {
            fs::create_dir(&dir)?;
        }
        Ok(dir)
    }
}
