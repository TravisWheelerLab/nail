use crate::align::bounded::structs::bound::CloudBoundGroup;
use crate::align::bounded::structs::row_bound_params::RowBoundParams;
use crate::structs::profile::constants::{
  SPECIAL_B, SPECIAL_C, SPECIAL_E, SPECIAL_J, SPECIAL_LOOP, SPECIAL_MOVE, SPECIAL_N,
};
use crate::structs::{DpMatrix, Profile};

pub fn posterior_bounded(
  profile: &Profile,
  forward_matrix: &mut DpMatrix,
  backward_matrix: &mut DpMatrix,
  posterior_matrix: &mut DpMatrix,
  params: &RowBoundParams,
) {
  let target_length = forward_matrix.target_length;
    let overall_score: f32 = forward_matrix.get_special(target_length, SPECIAL_C)
        + profile.special_transition_score(SPECIAL_C, SPECIAL_MOVE);
    let mut denominator: f32;

    posterior_matrix.set_special(0, SPECIAL_E, 0.0);
    posterior_matrix.set_special(0, SPECIAL_N, 0.0);
    posterior_matrix.set_special(0, SPECIAL_J, 0.0);
    posterior_matrix.set_special(0, SPECIAL_B, 0.0);
    posterior_matrix.set_special(0, SPECIAL_C, 0.0);

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

        posterior_matrix.set_special(i, SPECIAL_E, 0.0);
        posterior_matrix.set_special(
            i,
            SPECIAL_N,
            (forward_matrix.get_special(i - 1, SPECIAL_N)
                + backward_matrix.get_special(i, SPECIAL_N)
                + profile.special_transition_score(SPECIAL_N, SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(
            i,
            SPECIAL_J,
            (forward_matrix.get_special(i - 1, SPECIAL_J)
                + backward_matrix.get_special(i, SPECIAL_J)
                + profile.special_transition_score(SPECIAL_J, SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(i, SPECIAL_B, 0.0);

        posterior_matrix.set_special(
            i,
            SPECIAL_C,
            (forward_matrix.get_special(i - 1, SPECIAL_C)
                + backward_matrix.get_special(i, SPECIAL_C)
                + profile.special_transition_score(SPECIAL_C, SPECIAL_LOOP)
                - overall_score)
                .exp(),
        );

        denominator += posterior_matrix.get_special(i, SPECIAL_N);
        denominator += posterior_matrix.get_special(i, SPECIAL_J);
        denominator += posterior_matrix.get_special(i, SPECIAL_C);

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
            SPECIAL_N,
            posterior_matrix.get_special(i, SPECIAL_N) * denominator,
        );
        posterior_matrix.set_special(
            i,
            SPECIAL_J,
            posterior_matrix.get_special(i, SPECIAL_J) * denominator,
        );
        posterior_matrix.set_special(
            i,
            SPECIAL_C,
            posterior_matrix.get_special(i, SPECIAL_C) * denominator,
        );
    }
}
