use std::{
    cmp::Ordering::{Equal, Greater, Less},
    fmt::{Debug, Display},
    io::{BufWriter, Write},
    ops::{AddAssign, Div, DivAssign, MulAssign, SubAssign},
};

use lazy_static::lazy_static;

pub trait MaxAssign {
    fn max_assign(&mut self, other: Self);
}

pub trait MinAssign {
    fn min_assign(&mut self, other: Self);
}

impl<T: PartialOrd + Copy> MaxAssign for T {
    fn max_assign(&mut self, other: Self) {
        if *self < other {
            *self = other;
        }
    }
}

impl<T: PartialOrd + Copy> MinAssign for T {
    fn min_assign(&mut self, other: Self) {
        if *self > other {
            *self = other;
        }
    }
}

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

pub trait IterPrint {
    fn print_each(&self);
    fn write_display<W: Write>(&self, out: W) -> anyhow::Result<()>;
}

pub trait IterDebug {
    fn debug_each(&self);
    fn write_debug<W: Write>(&self, out: W) -> anyhow::Result<()>;
}

impl<I, T> IterPrint for I
where
    I: IntoIterator<Item = T> + ?Sized,
    for<'a> &'a I: IntoIterator<Item = &'a T>,
    T: Display,
{
    fn print_each(&self) {
        self.into_iter().for_each(|i| println!("{i}"));
    }

    fn write_display<W: Write>(&self, out: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(out);
        for (i, item) in self.into_iter().enumerate() {
            writeln!(out, "{i}: {item}")?;
        }
        Ok(())
    }
}

impl<I, T> IterDebug for I
where
    I: IntoIterator<Item = T> + ?Sized,
    for<'a> &'a I: IntoIterator<Item = &'a T>,
    T: Debug,
{
    fn debug_each(&self) {
        self.into_iter().for_each(|i| println!("{i:?}"));
    }

    fn write_debug<W: Write>(&self, out: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(out);
        for (i, item) in self.into_iter().enumerate() {
            writeln!(out, "{i}: {item:?}")?;
        }
        Ok(())
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
    fn mean(&self) -> Option<T>;
    fn argmax(&self) -> Option<usize>;
}

impl<T, U> VecMath<T> for U
where
    T: Float,
    U: AsRef<[T]>,
{
    fn mean(&self) -> Option<T> {
        let mut sum = T::from_usize(0);
        self.as_ref().iter().for_each(|&item| sum += item);

        Some(sum / T::from_usize(self.as_ref().len()))
    }

    fn argmax(&self) -> Option<usize> {
        let mut max = *self.as_ref().first()?;
        let mut argmax: usize = 0;

        for (idx, &item) in self.as_ref().iter().enumerate().skip(1) {
            if item > max {
                max = item;
                argmax = idx;
            }
        }

        Some(argmax)
    }
}

pub trait VecMathMut<T>
where
    T: Float,
{
    fn normalize(&mut self);
    fn add(&mut self, other: &[T]);
    fn sub(&mut self, other: &[T]);
    fn scale(&mut self, factor: T);
    fn saturate_lower(&mut self, min: T);
}

impl<T, U> VecMathMut<T> for U
where
    T: Float,
    U: AsMut<[T]>,
{
    fn normalize(&mut self) {
        let mut sum = T::from_usize(0);
        self.as_mut().iter().for_each(|&item| sum += item);
        self.as_mut().iter_mut().for_each(|item| *item /= sum);
    }

    fn add(&mut self, other: &[T]) {
        self.as_mut()
            .iter_mut()
            .zip(other)
            .for_each(|(a, &b)| *a += b);
    }

    fn sub(&mut self, other: &[T]) {
        self.as_mut()
            .iter_mut()
            .zip(other)
            .for_each(|(a, &b)| *a -= b);
    }

    fn scale(&mut self, factor: T) {
        self.as_mut().iter_mut().for_each(|item| *item *= factor);
    }

    fn saturate_lower(&mut self, min: T) {
        self.as_mut().iter_mut().for_each(|item| {
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
        .mean()
        .unwrap()
}

pub fn relative_entropy(a: &[f32], b: &[f32]) -> f32 {
    a.iter()
        .zip(b)
        .map(|(p_a, p_b)| p_a * (p_a / p_b).log2())
        .sum()
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
