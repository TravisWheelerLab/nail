use anyhow::Result;
use nale::alignment::{backward, forward, optimal_accuracy, posterior};
use nale::structs::{DpMatrix, Hmm, Profile, Sequence};
use nale::util::Dump;
use std::fs::File;
use std::io::BufWriter;

fn main() -> Result<()> {
    // let path = format!("{}/{}", env!("CARGO_MANIFEST_DIR"), "Pfam-A.hmm");
    let hmm_path = format!(
        "{}/{}",
        env!("CARGO_MANIFEST_DIR"),
        "resources/Alpha-amylase.hmm"
    );
    let seq_path = format!(
        "{}/{}",
        env!("CARGO_MANIFEST_DIR"),
        "resources/Alpha-amylase.c.fa"
    );

    let hmm_list = Hmm::from_p7hmm_file(hmm_path)?;
    let seq_list = Sequence::from_fasta(seq_path)?;

    // TODO: fix parser so it doesn't give an uninitialized model at the end of the list
    // println!("hmm count: {}", hmm_list.len() - 1);
    // println!("seq count: {}", seq_list.len());

    for hmm_idx in 0..hmm_list.len() - 1 {
        let mut profile = Profile::new(&hmm_list[hmm_idx]);
        for seq_idx in 0..seq_list.len() {
            let target = &seq_list[seq_idx];
            profile.configure_for_length(target.length);

            // TODO: reuse matrices
            let mut forward_matrix = DpMatrix::new(profile.length, target.length);
            forward(&profile, target, &mut forward_matrix);

            let mut forward_out = BufWriter::new(File::create("./nale-dump/forward.mtx")?);
            forward_matrix.dump(&mut forward_out)?;

            let mut backward_matrix = DpMatrix::new(profile.length, target.length);
            backward(&profile, target, &mut backward_matrix);

            let mut backward_out = BufWriter::new(File::create("./nale-dump/backward.mtx")?);
            backward_matrix.dump(&mut backward_out)?;

            let mut posterior_matrix = DpMatrix::new(profile.length, target.length);
            posterior(
                &profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
            );

            let mut posterior_out = BufWriter::new(File::create("./nale-dump/posterior.mtx")?);
            posterior_matrix.dump(&mut posterior_out)?;

            let mut optimal_matrix = DpMatrix::new(profile.length, target.length);
            optimal_accuracy(&profile, &posterior_matrix, &mut optimal_matrix);

            let mut optimal_out = BufWriter::new(File::create("./nale-dump/optimal.mtx")?);
            optimal_matrix.dump(&mut optimal_out)?;
            
            
        }
    }

    Ok(())
}
