use anyhow::Context;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

pub trait PathBufExt {
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>>;
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
}
