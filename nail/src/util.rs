use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::Context;
use image::{Rgb, RgbImage};
use libnail::align::structs::{AntiDiagonal, AntiDiagonalBounds, RowBounds};
use thiserror::Error;

#[derive(Default, Debug, Clone)]
pub enum FileFormat {
    Fasta,
    Stockholm,
    Hmm,
    #[default]
    Unset,
}

impl Display for FileFormat {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::Fasta => write!(f, "Fasta"),
            FileFormat::Stockholm => write!(f, "Stockholm"),
            FileFormat::Hmm => write!(f, "HMM"),
            FileFormat::Unset => write!(f, "Unset"),
        }
    }
}

#[derive(Error, Debug)]
#[error("can't guess file format of: {path}")]
pub struct UnrecognizedFileFormatError {
    path: String,
}

pub fn guess_query_format_from_query_file(
    query_path: &impl AsRef<Path>,
) -> anyhow::Result<FileFormat> {
    let file = File::open(query_path).context(format!(
        "failed to open query file: {}",
        query_path.as_ref().to_string_lossy()
    ))?;

    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    if &first_line[0..1] == ">" {
        Ok(FileFormat::Fasta)
    } else if &first_line[0..11] == "# STOCKHOLM" {
        Ok(FileFormat::Stockholm)
    } else if &first_line[0..5] == "HMMER" {
        Ok(FileFormat::Hmm)
    } else {
        Err(UnrecognizedFileFormatError {
            path: query_path.as_ref().to_string_lossy().to_string(),
        }
        .into())
    }
}

#[derive(Error, Debug)]
#[error("command exited without success")]
struct CommandExitStatusError;

/// An extension trait that is intended to add a run method to the std::process::Command struct.
pub trait CommandExt {
    fn run(&mut self) -> anyhow::Result<()>;
}

impl CommandExt for Command {
    fn run(&mut self) -> anyhow::Result<()> {
        let output = self.output().context("failed to run command")?;

        match output.status.success() {
            true => Ok(()),
            false => {
                let stdout = std::str::from_utf8(&output.stdout)
                    .context("failed to convert sdtout to UTF8")?;
                let stderr = std::str::from_utf8(&output.stderr)
                    .context("failed to convert sdterr to UTF8")?;

                println!("command:\n{self:?}\n");
                println!("stdout:\n{stdout}\n");
                println!("stderr:\n{stderr}\n");
                Err(CommandExitStatusError.into())
            }
        }
    }
}

pub trait PathBufExt {
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>>;
    fn remove(&self) -> anyhow::Result<()>;
}

impl PathBufExt for PathBuf {
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>> {
        let mut file_options = File::options();

        if allow_overwrite {
            file_options.write(true).truncate(true).create(true);
        } else {
            file_options.write(true).create_new(true);
        };

        let file = file_options
            .open(self)
            .context(format!("failed to create file: {}", self.to_string_lossy()))?;

        Ok(BufWriter::new(file))
    }

    fn remove(&self) -> anyhow::Result<()> {
        std::fs::remove_file(self)?;
        Ok(())
    }
}

pub fn cloud_image(
    bounds_a: &AntiDiagonalBounds,
    bounds_b: &AntiDiagonalBounds,
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let red = Rgb([255, 0, 0]);
    let green = Rgb([0, 0, 255]);
    let blue = Rgb([0, 255, 0]);

    let width = bounds_a.profile_length;
    let height = bounds_a.target_length;

    let mut img = RgbImage::new(width as u32 + 1, height as u32 + 1);

    let mut draw = |bound: &AntiDiagonal, color: Rgb<u8>| {
        bound.cell_zip().for_each(|(y, x)| {
            img.put_pixel(x as u32, y as u32, color);
        });
    };

    bounds_a.bounds().iter().for_each(|b| draw(b, red));
    bounds_b.bounds().iter().for_each(|b| draw(b, green));

    if let Some(interval) = AntiDiagonalBounds::intersection(bounds_a, bounds_b) {
        (interval.start..=interval.end)
            .map(|idx| (bounds_a.get(idx), bounds_b.get(idx)))
            .for_each(|(b_a, b_b)| {
                draw(b_a, blue);
                draw(b_b, blue);
            });
    }

    img.save(path).context("failed to write image")
}

#[allow(dead_code)]
pub trait Image {
    fn image(&self, path: impl AsRef<Path>) -> anyhow::Result<()>;
}

impl Image for AntiDiagonalBounds {
    fn image(&self, path: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut width = self.profile_length;
        let mut height = self.target_length;

        self.bounds().iter().for_each(|b| {
            width = width.max(b.left_profile_idx);
            width = width.max(b.right_profile_idx);
            height = height.max(b.left_target_idx);
            height = height.max(b.right_target_idx);
        });

        let mut img = RgbImage::new(width as u32 + 1, height as u32 + 1);

        (0..=self.profile_length)
            .for_each(|x| img.put_pixel(x as u32, self.target_length as u32, Rgb([0, 255, 0])));

        (0..=self.target_length)
            .for_each(|y| img.put_pixel(self.profile_length as u32, y as u32, Rgb([0, 255, 0])));

        for (idx, bound) in self.bounds.iter().enumerate() {
            let color = if idx % 2 == 0 {
                Rgb([255, 0, 0])
            } else {
                Rgb([0, 0, 255])
            };

            bound.cell_zip().for_each(|(y, x)| {
                img.put_pixel(x as u32, y as u32, color);
            });
        }
        img.save(path).context("failed to write image")
    }
}

impl Image for RowBounds {
    fn image(&self, path: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut width = self.profile_length;
        let height = self.target_length;

        (self.target_start..=self.target_end)
            .map(|idx| (self.left_row_bounds[idx], self.right_row_bounds[idx]))
            .for_each(|(left, right)| {
                if left != usize::MAX {
                    width = width.max(left);
                }
                width = width.max(right);
            });

        let mut img = RgbImage::new(width as u32 + 1, height as u32 + 1);

        (0..=self.profile_length)
            .for_each(|x| img.put_pixel(x as u32, self.target_length as u32, Rgb([0, 255, 0])));

        (0..=self.target_length)
            .for_each(|y| img.put_pixel(self.profile_length as u32, y as u32, Rgb([0, 255, 0])));

        (self.target_start..=self.target_end)
            .map(|idx| (self.left_row_bounds[idx], self.right_row_bounds[idx]))
            .enumerate()
            .for_each(|(y, (left, right))| {
                if left != usize::MAX {
                    (left..=right).for_each(|x| {
                        img.put_pixel(x as u32, y as u32, Rgb([255, 0, 0]));
                    })
                }
            });

        img.save(path).context("failed to write image")
    }
}

pub fn check_hmmer_installed() -> anyhow::Result<()> {
    Command::new("hmmbuild")
        .arg("-h")
        .run()
        .context("hmmbuild does not appear to be in the system path")
}

pub fn check_mmseqs_installed() -> anyhow::Result<()> {
    Command::new("mmseqs")
        .arg("-h")
        .run()
        .context("mmseqs2 does not appear to be in the system path")
}

pub fn set_threads(num_threads: usize) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("failed to build rayon global threadpool")
}
