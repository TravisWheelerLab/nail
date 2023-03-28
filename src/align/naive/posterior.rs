use crate::structs::dp_matrix::DpMatrix;
use crate::structs::Profile;
use crate::timing::time;

#[funci::timed(timer = time)]
pub fn posterior(
    profile: &Profile,
    forward_matrix: &impl DpMatrix,
    backward_matrix: &impl DpMatrix,
    posterior_matrix: &mut impl DpMatrix,
) {
    let target_length = forward_matrix.target_length();
    let overall_score: f32 = forward_matrix.get_special(target_length, Profile::SPECIAL_C)
        + profile.special_transition_score(Profile::SPECIAL_C, Profile::SPECIAL_MOVE);
    let mut denominator: f32;

    posterior_matrix.set_special(0, Profile::SPECIAL_E, 0.0);
    posterior_matrix.set_special(0, Profile::SPECIAL_N, 0.0);
    posterior_matrix.set_special(0, Profile::SPECIAL_J, 0.0);
    posterior_matrix.set_special(0, Profile::SPECIAL_B, 0.0);
    posterior_matrix.set_special(0, Profile::SPECIAL_C, 0.0);

    for k in 0..=profile.length {
        posterior_matrix.set_match(0, k, 0.0);
        posterior_matrix.set_insert(0, k, 0.0);
        posterior_matrix.set_delete(0, k, 0.0);
    }

    for i in 1..=target_length {
        denominator = 0.0;
        posterior_matrix.set_match(i, 0, 0.0);
        posterior_matrix.set_insert(i, 0, 0.0);
        posterior_matrix.set_delete(i, 0, 0.0);

        for k in 1..profile.length {
            posterior_matrix.set_match(
                i,
                k,
                (forward_matrix.get_match(i, k) + backward_matrix.get_match(i, k) - overall_score)
                    .exp(),
            );

            denominator += posterior_matrix.get_match(i, k);

            posterior_matrix.set_insert(
                i,
                k,
                (forward_matrix.get_insert(i, k) + backward_matrix.get_insert(i, k)
                    - overall_score)
                    .exp(),
            );
            denominator += posterior_matrix.get_insert(i, k);

            posterior_matrix.set_delete(i, k, 0.0);
        }

        posterior_matrix.set_match(
            i,
            profile.length,
            (forward_matrix.get_match(i, profile.length)
                + backward_matrix.get_match(i, profile.length)
                - overall_score)
                .exp(),
        );

        denominator += posterior_matrix.get_match(i, profile.length);
        posterior_matrix.set_insert(i, profile.length, 0.0);
        posterior_matrix.set_delete(i, profile.length, 0.0);

        posterior_matrix.set_special(i, Profile::SPECIAL_E, 0.0);
        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_N,
            (forward_matrix.get_special(i - 1, Profile::SPECIAL_N)
                + backward_matrix.get_special(i, Profile::SPECIAL_N)
                + profile.special_transition_score(Profile::SPECIAL_N, Profile::SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_J,
            (forward_matrix.get_special(i - 1, Profile::SPECIAL_J)
                + backward_matrix.get_special(i, Profile::SPECIAL_J)
                + profile.special_transition_score(Profile::SPECIAL_J, Profile::SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(i, Profile::SPECIAL_B, 0.0);

        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_C,
            (forward_matrix.get_special(i - 1, Profile::SPECIAL_C)
                + backward_matrix.get_special(i, Profile::SPECIAL_C)
                + profile.special_transition_score(Profile::SPECIAL_C, Profile::SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        denominator += posterior_matrix.get_special(i, Profile::SPECIAL_N);
        denominator += posterior_matrix.get_special(i, Profile::SPECIAL_J);
        denominator += posterior_matrix.get_special(i, Profile::SPECIAL_C);

        denominator = 1.0 / denominator;

        for k in 1..profile.length {
            posterior_matrix.set_match(i, k, posterior_matrix.get_match(i, k) * denominator);
            posterior_matrix.set_insert(i, k, posterior_matrix.get_insert(i, k) * denominator);
        }
        posterior_matrix.set_match(
            i,
            profile.length,
            posterior_matrix.get_match(i, profile.length) * denominator,
        );
        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_N,
            posterior_matrix.get_special(i, Profile::SPECIAL_N) * denominator,
        );
        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_J,
            posterior_matrix.get_special(i, Profile::SPECIAL_J) * denominator,
        );
        posterior_matrix.set_special(
            i,
            Profile::SPECIAL_C,
            posterior_matrix.get_special(i, Profile::SPECIAL_C) * denominator,
        );
    }
}
