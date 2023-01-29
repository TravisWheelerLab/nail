use crate::align::bounded::cloud_search;
use crate::align::bounded::structs::CloudMatrix;
use crate::structs::{Profile, Sequence};
use anyhow::Result;

pub fn pipeline_bounded(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    for profile in profiles.iter_mut() {
        for target in targets.iter() {
            profile.configure_for_length(target.length);

            let mut cloud_matrix = CloudMatrix::new(profile.length, target.length);
            cloud_search(profile, target, &mut cloud_matrix)?;
        }
    }

    Ok(())
}
