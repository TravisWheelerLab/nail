use crate::align::bounded::cloud_search;
use crate::align::bounded::structs::{CloudMatrix, CloudSearchParams};
use crate::structs::{Profile, Sequence};
use anyhow::Result;

pub fn pipeline_bounded(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    let mut cloud_search_params = CloudSearchParams::default();

    for profile in profiles.iter_mut() {
        for target in targets.iter() {
            profile.configure_for_length(target.length);

            cloud_search_params.target_start = 1;
            cloud_search_params.profile_start = 1;
            cloud_search_params.target_end = target.length;
            cloud_search_params.profile_end = profile.length;

            let mut cloud_matrix = CloudMatrix::new(profile.length, target.length);
            cloud_search(profile, target, &mut cloud_matrix, &cloud_search_params)?;
        }
    }

    Ok(())
}
