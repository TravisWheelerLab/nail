use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::SeedableRng;
use rand_pcg::Pcg64;

use libnail::{
    align::{
        cloud_search_bwd, cloud_search_fwd,
        structs::{AdMatrixLinear, Cloud, Seed},
        CloudSearchParams,
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

fn bench_cloud_search(c: &mut Criterion) {
    let params = CloudSearchParams::default();

    // Several (profile_len, seq_len) sizes to cover realistic range
    let cases: &[(usize, usize)] = &[
        (100, 200),
        (200, 500),
        (300, 1000),
        (500, 2000),
    ];

    let mut fwd_group = c.benchmark_group("cloud_search_fwd");
    for &(prf_len, seq_len) in cases {
        let (prf, seq) = make_fixture(prf_len, seq_len, 3, 5);
        let seed = Seed {
            prf: "bench".to_string(),
            seq: "bench".to_string(),
            prf_start: prf_len / 4,
            prf_end: 3 * prf_len / 4,
            seq_start: seq_len / 4,
            seq_end: 3 * seq_len / 4,
            score: 100.0,
            e_value: 1e-10,
            cigar: vec![],
            trace_bounds: Default::default(),
        };

        fwd_group.bench_with_input(
            BenchmarkId::new("prf_len", format!("{prf_len}x{seq_len}")),
            &(prf_len, seq_len),
            |b, _| {
                let mut mx = AdMatrixLinear::new(seq.length);
                let mut cloud = Cloud::new(seq.length, prf.length);
                b.iter(|| {
                    mx.reuse(seq.length);
                    cloud.reuse(seq.length, prf.length);
                    cloud_search_fwd(&prf, &seq, &seed, &mut mx, &params, &mut cloud)
                });
            },
        );
    }
    fwd_group.finish();

    let mut bwd_group = c.benchmark_group("cloud_search_bwd");
    for &(prf_len, seq_len) in cases {
        let (prf, seq) = make_fixture(prf_len, seq_len, 3, 5);
        let seed = Seed {
            prf: "bench".to_string(),
            seq: "bench".to_string(),
            prf_start: prf_len / 4,
            prf_end: 3 * prf_len / 4,
            seq_start: seq_len / 4,
            seq_end: 3 * seq_len / 4,
            score: 100.0,
            e_value: 1e-10,
            cigar: vec![],
            trace_bounds: Default::default(),
        };

        bwd_group.bench_with_input(
            BenchmarkId::new("prf_len", format!("{prf_len}x{seq_len}")),
            &(prf_len, seq_len),
            |b, _| {
                let mut mx = AdMatrixLinear::new(seq.length);
                let mut cloud = Cloud::new(seq.length, prf.length);
                b.iter(|| {
                    cloud_search_bwd(&prf, &seq, &seed, &mut mx, &params, &mut cloud)
                });
            },
        );
    }
    bwd_group.finish();
}

criterion_group!(benches, bench_cloud_search);
criterion_main!(benches);
