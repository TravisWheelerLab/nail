use std::fmt::Display;

/// A single CIGAR operation from the mmseqs alignment trace.
/// Convention follows mmseqs/SAM: operations are relative to the query (profile).
#[derive(Clone, Debug, PartialEq)]
pub enum CigarOp {
    /// Match/mismatch: both profile and sequence advance
    Match(usize),
    /// Insertion in query: profile advances, sequence does not (gap in sequence)
    Insertion(usize),
    /// Deletion in query: sequence advances, profile does not (gap in profile)
    Deletion(usize),
}

/// Parse a CIGAR string into a vector of CigarOps.
pub fn parse_cigar(cigar: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num = 0usize;
    for c in cigar.bytes() {
        match c {
            b'0'..=b'9' => num = num * 10 + (c - b'0') as usize,
            b'M' => { ops.push(CigarOp::Match(num)); num = 0; }
            b'I' => { ops.push(CigarOp::Insertion(num)); num = 0; }
            b'D' => { ops.push(CigarOp::Deletion(num)); num = 0; }
            _ => { num = 0; }
        }
    }
    ops
}

/// Precomputed per-AD bounds from the mmseqs alignment trace.
///
/// For each anti-diagonal that the trace crosses, stores the minimum
/// profile range (prf_lo..prf_hi) that must be included in the cloud
/// bound, expanded by ±BAND cells around the trace position.
///
/// Built once from the CIGAR at seed parse time. Indexed by AD number.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TraceBounds {
    /// prf_lo[ad] = minimum profile index that must be in the bound (0 = no constraint)
    pub prf_lo: Vec<usize>,
    /// prf_hi[ad] = maximum profile index that must be in the bound (0 = no constraint)
    pub prf_hi: Vec<usize>,
}

impl TraceBounds {
    const BAND: usize = 3;

    pub fn from_cigar(
        cigar: &[CigarOp],
        prf_start: usize,
        seq_start: usize,
        prf_len: usize,
        seq_len: usize,
    ) -> Self {
        if cigar.is_empty() {
            return Self::default();
        }

        let max_ad = prf_len + seq_len;
        let mut prf_lo = vec![0usize; max_ad + 1];
        let mut prf_hi = vec![0usize; max_ad + 1];

        let mut prf = prf_start;
        let mut seq = seq_start;

        // Clamp bounds to leave room for DP lookups at +1 offsets
        let prf_max = prf_len.saturating_sub(1);
        let seq_max = seq_len.saturating_sub(1);

        let mut set_bounds = |prf: usize, seq: usize, prf_lo: &mut Vec<usize>, prf_hi: &mut Vec<usize>| {
            let ad = prf + seq;
            if ad <= max_ad && ad < prf_lo.len() {
                // Clamp so that the bound doesn't push seq or prf past safe limits
                let lo = prf.saturating_sub(Self::BAND).max(1);
                let hi = (prf + Self::BAND).min(prf_max);
                // Also ensure seq stays in range: seq = ad - prf, so prf_lo >= ad - seq_max
                let lo = lo.max(ad.saturating_sub(seq_max));
                prf_lo[ad] = lo;
                prf_hi[ad] = hi;
            }
        };

        set_bounds(prf, seq, &mut prf_lo, &mut prf_hi);

        for op in cigar {
            let (count, dp, ds) = match op {
                CigarOp::Match(n) => (*n, 1usize, 1usize),
                CigarOp::Insertion(n) => (*n, 1, 0),
                CigarOp::Deletion(n) => (*n, 0, 1),
            };
            for _ in 0..count {
                prf += dp;
                seq += ds;
                set_bounds(prf, seq, &mut prf_lo, &mut prf_hi);
            }
        }

        // Fill gaps: Match ops advance AD by 2, leaving odd ADs empty.
        // Interpolate from neighbors so there are no holes in coverage.
        for ad in 1..max_ad {
            if prf_hi[ad] == 0 {
                let has_prev = prf_hi[ad.saturating_sub(1)] > 0;
                let has_next = ad + 1 <= max_ad && prf_hi[ad + 1] > 0;
                if has_prev && has_next {
                    prf_lo[ad] = prf_lo[ad - 1].min(prf_lo[ad + 1]);
                    prf_hi[ad] = prf_hi[ad - 1].max(prf_hi[ad + 1]);
                }
            }
        }

        Self { prf_lo, prf_hi }
    }

    /// Returns true if the trace covers this anti-diagonal.
    #[inline]
    pub fn has(&self, ad: usize) -> bool {
        ad < self.prf_hi.len() && self.prf_hi[ad] > 0
    }

    /// Get the required profile range for this AD, or None.
    #[inline]
    pub fn get(&self, ad: usize) -> Option<(usize, usize)> {
        if self.has(ad) {
            Some((self.prf_lo[ad], self.prf_hi[ad]))
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Seed {
    pub prf: String,
    pub seq: String,
    pub seq_start: usize,
    pub seq_end: usize,
    pub prf_start: usize,
    pub prf_end: usize,
    pub score: f32,
    pub e_value: f64,
    pub cigar: Vec<CigarOp>,
    pub trace_bounds: TraceBounds,
}

impl Display for Seed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}:{}",
            self.seq_start, self.seq_end, self.prf_start, self.prf_end, self.score
        )
    }
}
