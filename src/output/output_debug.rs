use anyhow::Context;
use std::fs::{create_dir_all, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

pub fn get_profile_target_output_dir_path(
    root_path: &Path,
    profile_name: &str,
    target_name: &str,
) -> anyhow::Result<PathBuf> {
    let path = root_path.join(profile_name).join(target_name);
    create_dir_all(&path)?;
    Ok(path)
}

pub fn set_file_name_and_get_buf_writer(
    path: &mut PathBuf,
    file_name: &str,
    allow_overwrite: bool,
) -> anyhow::Result<BufWriter<File>> {
    if path.is_file() {
        path.set_file_name(file_name);
    } else {
        path.push(file_name);
    }

    let mut file_options = File::options();

    if allow_overwrite {
        file_options.write(true).truncate(true).create(true);
    } else {
        file_options.write(true).create_new(true);
    };

    let file = file_options
        .open(&path)
        .context(format!("failed to create file: {}", path.to_string_lossy()))?;

    Ok(BufWriter::new(file))
}
