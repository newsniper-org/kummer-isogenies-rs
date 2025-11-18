use criterion::{criterion_group, criterion_main, Criterion};
use cubemoma::FpComplex;
use kummer_isogenies_rs::{biguint_to_fp, field_types::{Fp, Fp2}, fp_to_biguint, kuhash::KummerHashSecp256, params::PARAMS};
use num_bigint::BigUint;
use rand::{CryptoRng, Rng};
use std::{fmt::format, hint::black_box};

fn random_msg<R: CryptoRng + Sized>(rng: &mut R) -> [Fp; 3] {
    let msgs: [Fp; 3] = [0, 1, 2].map(|_| {
        let mut rnd: [u64; 4] = [0u64; 4];
        loop {
            for limb in rnd.iter_mut() {
                *limb = rng.random();
            }
            let interm = Fp::new(rnd);
            if interm == Fp::zero() {
                continue;
            } else if interm < biguint_to_fp(BigUint::from(3u64).pow(PARAMS.e - 1)) {
                break;
            } else {
                continue;
            }
        }
        let (q, r) = (PARAMS.l as usize / 64, PARAMS.l as usize % 64);
        rnd[q] += 1u64 << r;
        Fp::new(rnd)
    });
    msgs
}


fn benchmark_hash() -> [Fp2; 4] {
    let hash = KummerHashSecp256::new();
    let [m1, m2, m3] = random_msg(&mut rand::rng());

    let h = hash.hash_optimal(&[m1, m2, m3]);
    h
}


fn bench_cpu_kuhash(c: &mut Criterion) {
    c.bench_function("cpu_kuhash", |bencher| {
        bencher.iter(|| black_box(benchmark_hash()));
    });
}

criterion_group!(benches, bench_cpu_kuhash);
criterion_main!(benches);