use crate::structs::Profile;

use super::DpMatrix;

#[derive(Default, Clone)]
pub struct DpMatrixFlat {
    pub target_length: usize,
    pub profile_length: usize,
    /// The DP matrix core model data cells as a flat vector.
    ///
    /// It's stored in the following pattern:
    ///     [
    ///
    ///         m_(0, 0), i_(0, 0), d_(0, 0),
    ///         m_(0, 1), i_(0, 1), d_(0, 1),
    ///         ...
    ///         m_(0, P), i_(0, P), d_(0, P),
    ///         ...
    ///         m_(T, 0), i_(T, 0), d_(T, 0),
    ///         ...
    ///         m_(T, P), i_(T, P), d_(T, P)
    ///
    ///     ]
    ///
    /// where:
    ///
    ///     T:        <target_length>
    ///     S:        <profile_length>
    ///     m_(i, j): the match score at cell (i, j)
    ///     i_(i, j): the insert score at cell (i, j)
    ///     d_(i, j): the delete score at cell (i, j)
    ///
    pub core_data: Vec<f32>,
    /// The DP matrix special state data cells as a flat vector.
    ///
    /// It's stored in the following pattern:
    ///     [
    ///
    ///         N_0, B_0, E_0, C_0, J_0,
    ///         N_1, B_1, E_1, C_1, J_1,
    ///         ...
    ///         N_T, B_T, E_T, C_T, J_T,
    ///
    ///     ]
    ///
    /// where:
    ///
    ///     T:        <target_length>
    ///     N_i: the N state score at target position i
    ///     B_i: the B state score at target position i
    ///     E_i: the E state score at target position i
    ///     C_i: the C state score at target position i
    ///     J_i: the J state score at target position i
    ///
    pub special_data: Vec<f32>,
}

impl DpMatrixFlat {
    pub fn new(target_length: usize, profile_length: usize) -> Self {
        let core_length = 3 * (target_length + 1) * (profile_length + 1);
        let special_length = 5 * (target_length + 1);
        DpMatrixFlat {
            target_length,
            profile_length,
            core_data: vec![-f32::INFINITY; core_length],
            special_data: vec![-f32::INFINITY; special_length],
        }
    }

    pub fn resize(&mut self, new_target_length: usize, new_profile_length: usize) {
        let new_size = new_target_length * new_profile_length;
        if new_size > self.core_data.len() {
            self.core_data.resize(new_size, -f32::INFINITY);
        }
        self.target_length = new_target_length;
        self.profile_length = new_profile_length;
    }

    pub fn reset(&mut self) {
        for core_idx in 0..(3 * (self.target_length + 1) * (self.profile_length + 1)) {
            self.core_data[core_idx] = -f32::INFINITY;
        }

        for special_idx in 0..((self.target_length + 1) * 5) {
            self.special_data[special_idx] = -f32::INFINITY;
        }
    }

    pub fn reuse(&mut self, new_target_length: usize, new_profile_length: usize) {
        // TODO: need logic for resizing self.special_data
        let new_core_length = 3 * (new_target_length + 1) * (new_profile_length + 1);
        let new_special_length = 5 * (new_target_length + 1);

        if new_core_length > self.core_data.len() {
            panic!("tried to resize DpMatrix: core");
        } else if new_special_length > self.special_data.len() {
            panic!("tried to resize DpMatrix: special");
        }

        self.target_length = new_target_length;
        self.profile_length = new_profile_length;

        self.reset();
    }
}

impl DpMatrix for DpMatrixFlat {
    fn target_length(&self) -> usize {
        self.target_length
    }

    fn profile_length(&self) -> usize {
        self.profile_length
    }

    #[inline]
    fn get_match(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (3 * profile_idx)]
    }

    #[inline]
    fn set_match(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + (3 * profile_idx)] = value;
    }

    #[inline]
    fn get_insert(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 1)]
    }

    #[inline]
    fn set_insert(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 1)] =
            value;
    }

    #[inline]
    fn get_delete(&self, target_idx: usize, profile_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 2)]
    }

    #[inline]
    fn set_delete(&mut self, target_idx: usize, profile_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(profile_idx <= self.profile_length);
        self.core_data[target_idx * (3 * (self.profile_length + 1)) + ((3 * profile_idx) + 2)] =
            value;
    }

    #[inline]
    fn get_special(&self, target_idx: usize, special_idx: usize) -> f32 {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx]
    }

    #[inline]
    fn set_special(&mut self, target_idx: usize, special_idx: usize, value: f32) {
        debug_assert!(target_idx <= self.target_length);
        debug_assert!(special_idx < Profile::NUM_SPECIAL_STATES);
        self.special_data[target_idx * 5 + special_idx] = value;
    }
}
