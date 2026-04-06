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

/// A cursor that walks along a CIGAR trace one anti-diagonal at a time.
/// Each CIGAR step advances the AD by 1. Call `advance_forward()` or
/// `advance_backward()` to step through the trace.
#[derive(Clone, Debug)]
pub struct TraceCursor {
    cigar: Vec<CigarOp>,
    op_idx: usize,
    op_offset: usize,
    pub prf_idx: usize,
    pub seq_idx: usize,
    pub active: bool,
}

impl TraceCursor {
    pub fn new_forward(cigar: &[CigarOp], prf_start: usize, seq_start: usize) -> Self {
        Self {
            cigar: cigar.to_vec(),
            op_idx: 0,
            op_offset: 0,
            prf_idx: prf_start,
            seq_idx: seq_start,
            active: !cigar.is_empty(),
        }
    }

    pub fn new_backward(cigar: &[CigarOp], prf_end: usize, seq_end: usize) -> Self {
        if cigar.is_empty() {
            return Self {
                cigar: Vec::new(), op_idx: 0, op_offset: 0,
                prf_idx: prf_end, seq_idx: seq_end, active: false,
            };
        }
        let last_op = cigar.len() - 1;
        let last_offset = match cigar[last_op] {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n) => n.saturating_sub(1),
        };
        Self {
            cigar: cigar.to_vec(),
            op_idx: last_op,
            op_offset: last_offset,
            prf_idx: prf_end,
            seq_idx: seq_end,
            active: true,
        }
    }

    pub fn ad(&self) -> usize {
        self.prf_idx + self.seq_idx
    }

    pub fn advance_forward(&mut self) -> bool {
        if !self.active { return false; }
        let op_len = match self.cigar[self.op_idx] {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n) => n,
        };
        self.op_offset += 1;
        if self.op_offset >= op_len {
            self.op_idx += 1;
            self.op_offset = 0;
            if self.op_idx >= self.cigar.len() {
                self.active = false;
                return false;
            }
        }
        match self.cigar[self.op_idx] {
            CigarOp::Match(_) => { self.prf_idx += 1; self.seq_idx += 1; }
            CigarOp::Insertion(_) => { self.prf_idx += 1; }
            CigarOp::Deletion(_) => { self.seq_idx += 1; }
        }
        true
    }

    pub fn advance_backward(&mut self) -> bool {
        if !self.active { return false; }
        match self.cigar[self.op_idx] {
            CigarOp::Match(_) => { self.prf_idx = self.prf_idx.saturating_sub(1); self.seq_idx = self.seq_idx.saturating_sub(1); }
            CigarOp::Insertion(_) => { self.prf_idx = self.prf_idx.saturating_sub(1); }
            CigarOp::Deletion(_) => { self.seq_idx = self.seq_idx.saturating_sub(1); }
        }
        if self.op_offset > 0 {
            self.op_offset -= 1;
        } else if self.op_idx > 0 {
            self.op_idx -= 1;
            let prev_len = match self.cigar[self.op_idx] {
                CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n) => n,
            };
            self.op_offset = prev_len.saturating_sub(1);
        } else {
            self.active = false;
        }
        true
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
