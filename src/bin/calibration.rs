//! Calibration of the replicable skew-symmetric share on simulated screens.
//!
//! Usage: calibration <screen1> <screen2> <umi_floor> [n_rep]   (TSV on standard output)
//!
//! A response matrix with a known skew-symmetric share is built from the observed data: its
//! symmetric part is the symmetric part of screen 2, and its skew-symmetric part is the
//! skew-symmetric part of screen 2 with gene labels permuted, rescaled to the target share.
//! Two screens are simulated by adding independent Gaussian noise with the standard errors
//! of the two real screens. The cross-screen estimator should return the target share, and
//! zero when the true matrix is symmetric, whereas the share computed within one screen is
//! inflated because independent noise divides equally between the two parts.
use anyhow::Result;
use perturbseq_reciprocity::{build, pearson, shared_genes, split, Screen};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

fn gauss(rng: &mut Xoshiro256PlusPlus) -> f64 {
    let (u, v): (f64, f64) = (rng.gen::<f64>().max(1e-300), rng.gen());
    (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}
fn dot(x: &[f64], y: &[f64]) -> f64 { x.iter().zip(y).map(|(a, b)| a * b).sum() }
fn upper(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1) / 2); for i in 0..n { for j in i + 1..n { v.push(x[i * n + j]); } } v }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let floor: f64 = a[3].parse()?;
    let nrep: usize = a.get(4).map(|s| s.parse().unwrap()).unwrap_or(20);
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let (ix1, ix2) = shared_genes(&d1, &d2, floor);
    let n = ix1.len();
    let m2 = build(&d2, &d2.mean, &ix2, 1.0).raw;
    let (s0, a_obs) = split(&m2, n);
    let se = |d: &Screen, ix: &[usize]| -> Vec<f64> { let mut v = vec![0f64; n * n]; for (x, &i) in ix.iter().enumerate() { for (y, &j) in ix.iter().enumerate() { v[x * n + y] = d.se(i, j); } } v };
    let (se1, se2) = (se(&d1, &ix1), se(&d2, &ix2));
    // skew part with permuted gene labels, unrelated to the symmetric part
    let mut perm: Vec<usize> = (0..n).collect();
    perm.shuffle(&mut Xoshiro256PlusPlus::seed_from_u64(7));
    let mut a0 = vec![0f64; n * n];
    for i in 0..n { for j in 0..n { a0[i * n + j] = a_obs[perm[i] * n + perm[j]]; } }
    let (ns, na) = (dot(&s0, &s0), dot(&a0, &a0));
    println!("umi_floor\tn_genes\ttrue_skew_share\tn_rep\tcross_screen_estimate_mean\tcross_screen_estimate_sd\tsingle_screen_share_mean\tr_skew_between_screens_mean");
    for &f in &[0.0f64, 0.1, 0.25, 0.5] {
        let c = if f == 0.0 { 0.0 } else { (f / (1.0 - f) * ns / na).sqrt() };
        let truth: Vec<f64> = s0.iter().zip(&a0).map(|(s, a)| s + c * a).collect();
        let res: Vec<(f64, f64, f64)> = (0..nrep).into_par_iter().map(|r| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(1000 + r as u64);
            let x1: Vec<f64> = truth.iter().zip(&se1).map(|(t, s)| t + s * gauss(&mut rng)).collect();
            let x2: Vec<f64> = truth.iter().zip(&se2).map(|(t, s)| t + s * gauss(&mut rng)).collect();
            let ((s1, a1), (s2, a2)) = (split(&x1, n), split(&x2, n));
            let (vs, va) = (dot(&s1, &s2), dot(&a1, &a2));
            (va / (vs + va), dot(&a1, &a1) / (dot(&s1, &s1) + dot(&a1, &a1)), pearson(&upper(&a1, n), &upper(&a2, n)))
        }).collect();
        let mean = |k: usize| res.iter().map(|r| [r.0, r.1, r.2][k]).sum::<f64>() / nrep as f64;
        let sd = (res.iter().map(|r| (r.0 - mean(0)).powi(2)).sum::<f64>() / (nrep as f64 - 1.0)).sqrt();
        println!("{floor}\t{n}\t{f}\t{nrep}\t{:.4}\t{:.4}\t{:.4}\t{:.4}", mean(0), sd, mean(1), mean(2));
    }
    Ok(())
}
