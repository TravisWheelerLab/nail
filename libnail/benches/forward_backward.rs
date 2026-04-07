use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::SeedableRng;
use rand_pcg::Pcg64;

use libnail::{
    align::{
        forward, backward, posterior,
        structs::{DpMatrixSparse, RowBounds},
    },
    structs::{Profile, Sequence},
};

fn make_fixture(prf_len: usize, seq_len: usize, prf_seed: u64, seq_seed: u64) -> (Profile, Sequence) {
    let mut rng = Pcg64::seed_from_u64(prf_seed);
    let prf_seq = Sequence::random_amino(prf_len, &mut rng);
    let mut prf = Profile::from_blosum_62_and_seq(&prf_seq).expect("valid profile");

    let mut rng2 = Pcg64::seed_from_u64(seq_seed);
    let seq = Sequence::random_amino(seq_len, &mut rng2);

    prf.configure_for_target_length(seq.length);
    (prf, seq)
}

fn full_bounds(prf: &Profile, seq: &Sequence) -> RowBounds {
    let mut bounds = RowBounds::new(seq.length);
    bounds.fill_rectangle(1, 1, seq.length, prf.length);
    bounds
}

fn bench_forward_backward(c: &mut Criterion) {
    let cases: &[(usize, usize)] = &[
        (100, 200),
        (200, 500),
        (300, 1000),
        (500, 2000),
    ];

    let mut fwd_group = c.benchmark_group("forward");
    for &(prf_len, seq_len) in cases {
        let (prf, seq) = make_fixture(prf_len, seq_len, 3, 5);
        let bounds = full_bounds(&prf, &seq);
        let mut mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);

        fwd_group.bench_with_input(
            BenchmarkId::new("prf_len", format!("{prf_len}x{seq_len}")),
            &(prf_len, seq_len),
            |b, _| {
                b.iter(|| {
                    mx.reuse(seq.length, prf.length, &bounds);
                    forward(&prf, &seq, &mut mx, &bounds)
                });
            },
        );
    }
    fwd_group.finish();

    let mut bwd_group = c.benchmark_group("backward");
    for &(prf_len, seq_len) in cases {
        let (prf, seq) = make_fixture(prf_len, seq_len, 3, 5);
        let bounds = full_bounds(&prf, &seq);
        let mut mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);

        bwd_group.bench_with_input(
            BenchmarkId::new("prf_len", format!("{prf_len}x{seq_len}")),
            &(prf_len, seq_len),
            |b, _| {
                b.iter(|| {
                    mx.reuse(seq.length, prf.length, &bounds);
                    backward(&prf, &seq, &mut mx, &bounds)
                });
            },
        );
    }
    bwd_group.finish();
}

fn bench_posterior(c: &mut Criterion) {
    let cases: &[(usize, usize)] = &[
        (100, 200),
        (200, 500),
        (300, 1000),
        (500, 2000),
    ];

    let mut group = c.benchmark_group("posterior");
    for &(prf_len, seq_len) in cases {
        let (prf, seq) = make_fixture(prf_len, seq_len, 3, 5);
        let bounds = full_bounds(&prf, &seq);
        let mut fwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let mut bwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let mut post_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);

        // Pre-fill fwd and bwd with real values
        forward(&prf, &seq, &mut fwd_mx, &bounds);
        backward(&prf, &seq, &mut bwd_mx, &bounds);

        group.bench_with_input(
            BenchmarkId::new("prf_len", format!("{prf_len}x{seq_len}")),
            &(prf_len, seq_len),
            |b, _| {
                b.iter(|| {
                    post_mx.reuse(seq.length, prf.length, &bounds);
                    posterior(&prf, &fwd_mx, &bwd_mx, &mut post_mx, &bounds)
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_forward_backward, bench_posterior);
criterion_main!(benches);
