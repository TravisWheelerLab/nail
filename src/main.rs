use anyhow::Result;
use clap::Parser;
use nale::align::{backward, forward, optimal_accuracy, posterior, traceback};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, DpMatrix, Profile, Sequence, Trace};
use std::fs::File;
use std::io::{stdout, BufWriter};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// query hmm
    query: String,
    /// target fasta
    target: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // println!("{:?}", args);

    let hmm_list = parse_hmms_from_p7hmm_file(args.query)?;
    let seq_list = Sequence::amino_from_fasta(&args.target)?;

    // println!("seq count: {}", seq_list.len());
    // println!("hmm count: {}", hmm_list.len());

    for hmm_idx in 0..hmm_list.len() {
        let mut profile = Profile::new(&hmm_list[hmm_idx]);
        for seq_idx in 0..seq_list.len() {
            let target = &seq_list[seq_idx];
            profile.configure_for_length(target.length);

            // TODO: reuse matrices
            let mut forward_matrix = DpMatrix::new(profile.length, target.length);
            forward(&profile, target, &mut forward_matrix);

            // let mut forward_out = BufWriter::new(File::create("./nale-dump/forward.mtx")?);
            // forward_matrix.dump(&mut forward_out)?;

            let mut backward_matrix = DpMatrix::new(profile.length, target.length);
            backward(&profile, target, &mut backward_matrix);

            // let mut backward_out = BufWriter::new(File::create("./nale-dump/backward.mtx")?);
            // backward_matrix.dump(&mut backward_out)?;

            let mut posterior_matrix = DpMatrix::new(profile.length, target.length);
            posterior(
                &profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
            );

            // let mut posterior_out = BufWriter::new(File::create("./nale-dump/posterior.mtx")?);
            // posterior_matrix.dump(&mut posterior_out)?;

            let mut optimal_matrix = DpMatrix::new(profile.length, target.length);
            optimal_accuracy(&profile, &posterior_matrix, &mut optimal_matrix);

            // let mut optimal_out = BufWriter::new(File::create("./nale-dump/optimal.mtx")?);
            // optimal_matrix.dump(&mut optimal_out)?;

            let mut trace = Trace::new(profile.length, target.length);
            traceback(&profile, &posterior_matrix, &optimal_matrix, &mut trace);

            let mut trace_out = BufWriter::new(File::create("./nale-dump/trace.dump")?);
            trace.dump(&mut trace_out, &profile, &target)?;

            let alignment = Alignment::new(&trace, &profile, &target);

            alignment.dump(&mut stdout())?;
        }
    }

    Ok(())
}
