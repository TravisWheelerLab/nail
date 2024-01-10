use crate::structs::Sequence;

pub enum SimpleTraceStep {
    Diagonal,
    Up,
    Left,
}

pub type SimpleTrace = Vec<SimpleTraceStep>;

#[allow(dead_code)]
pub fn print_trace(trace: &SimpleTrace) {
    for step in trace {
        match step {
            SimpleTraceStep::Diagonal => print!("D"),
            SimpleTraceStep::Up => print!("U"),
            SimpleTraceStep::Left => print!("L"),
        }
    }
    println!();
}

#[allow(dead_code)]
pub fn print_alignment(trace: &SimpleTrace, seq_1: &Sequence, seq_2: &Sequence) {
    let mut top = String::from("");
    let mut middle = String::from("");
    let mut bottom = String::from("");
    let mut seq_1_idx = 0;
    let mut seq_2_idx = 0;

    for step in trace {
        match step {
            SimpleTraceStep::Diagonal => {
                seq_1_idx += 1;
                seq_2_idx += 1;
                top.push(char::from(seq_1.utf8_bytes[seq_1_idx]));
                bottom.push(char::from(seq_2.utf8_bytes[seq_2_idx]));
                if seq_1.digital_bytes[seq_1_idx] == seq_2.digital_bytes[seq_2_idx] {
                    middle.push('|');
                } else {
                    middle.push(' ');
                }
            }
            SimpleTraceStep::Up => {
                seq_1_idx += 1;
                top.push(char::from(seq_1.utf8_bytes[seq_1_idx]));
                bottom.push('-');
                middle.push(' ');
            }
            SimpleTraceStep::Left => {
                seq_2_idx += 1;
                top.push('-');
                bottom.push(char::from(seq_2.utf8_bytes[seq_2_idx]));
                middle.push(' ');
            }
        }
    }
    println!("{}\n{}\n{}", top, middle, bottom)
}

const MATCH_SCORE: isize = 1;
const MISMATCH_SCORE: isize = -1;
const GAP_SCORE: isize = -2;

pub fn needleman_wunsch(seq_1: &Sequence, seq_2: &Sequence) -> SimpleTrace {
    let mut dp_matrix: Vec<Vec<isize>> = vec![vec![0; seq_2.length + 1]; seq_1.length + 1];

    for seq_2_idx in 0..=seq_2.length {
        dp_matrix[0][seq_2_idx] = (seq_2_idx as isize) * GAP_SCORE;
    }

    for seq_1_idx in 1..=seq_1.length {
        let seq_1_residue = seq_1.digital_bytes[seq_1_idx];
        dp_matrix[seq_1_idx][0] = (seq_1_idx as isize) * GAP_SCORE;

        for seq_2_idx in 1..=seq_2.length {
            let seq_2_residue = seq_2.digital_bytes[seq_2_idx];
            let match_score = if seq_1_residue == seq_2_residue {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            let diag_score = dp_matrix[seq_1_idx - 1][seq_2_idx - 1] + match_score;
            let up_score = dp_matrix[seq_1_idx - 1][seq_2_idx] + GAP_SCORE;
            let left_score = dp_matrix[seq_1_idx][seq_2_idx - 1] + GAP_SCORE;

            dp_matrix[seq_1_idx][seq_2_idx] = diag_score.max(up_score.max(left_score));
        }
    }

    let mut trace: SimpleTrace = vec![];
    let mut seq_1_idx = seq_1.length;
    let mut seq_2_idx = seq_2.length;

    while seq_1_idx > 0 && seq_2_idx > 0 {
        let seq_1_residue = seq_1.digital_bytes[seq_1_idx];
        let seq_2_residue = seq_2.digital_bytes[seq_2_idx];

        let match_score = if seq_1_residue == seq_2_residue {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        };

        let current_score = dp_matrix[seq_1_idx][seq_2_idx];

        let diag_target = current_score - match_score;
        let up_target = current_score - GAP_SCORE;
        let left_target = current_score - GAP_SCORE;

        let diag_score = dp_matrix[seq_1_idx - 1][seq_2_idx - 1];
        let up_score = dp_matrix[seq_1_idx - 1][seq_2_idx];
        let left_score = dp_matrix[seq_1_idx][seq_2_idx - 1];

        if diag_target == diag_score {
            seq_1_idx -= 1;
            seq_2_idx -= 1;
            trace.push(SimpleTraceStep::Diagonal);
        } else if up_target == up_score {
            seq_1_idx -= 1;
            trace.push(SimpleTraceStep::Up);
        } else if left_target == left_score {
            seq_2_idx -= 1;
            trace.push(SimpleTraceStep::Left);
        } else {
            panic!("simple traceback failed")
        }
    }

    while seq_1_idx > 0 {
        seq_1_idx -= 1;
        trace.push(SimpleTraceStep::Up);
    }

    while seq_2_idx > 0 {
        seq_2_idx -= 1;
        trace.push(SimpleTraceStep::Left);
    }

    trace.reverse();

    trace
}
