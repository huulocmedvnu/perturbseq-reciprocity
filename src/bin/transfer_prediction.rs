//! Prediction of the response matrix of one cell line from screens of other cell lines.
//!
//! Usage: transfer_prediction <umi_floor> <target> <reference1> <reference2> [...]   (TSV on standard output)
//!
//! Nothing but the control cells of the target line is used from the target screen. Every
//! prediction is compared with the response measured in the target line, for the whole
//! response and for its symmetric and skew-symmetric parts.
//!
//!   target_coexpression            c * rho of the target line
//!   reference_response             mean response matrix of the reference lines
//!   reference_adjusted             reference response + c * (rho of the target - rho of the reference)
//!   reference_skew_target_coexpr   c * rho of the target + skew-symmetric part of the reference response
//!   transfer_fitted                w1 * reference response + w2 * rho of the target + w3 * rho of the reference
//! The scale c is the least-squares slope of the reference response on the reference
//! co-expression. The weights w are fitted between pairs of reference lines.
//! target_half_vs_half is the agreement of two disjoint halves of the cells of the target screen.
use anyhow::{bail, Result};
use perturbseq_reciprocity::{pearson, read_f32, split, Screen};
use std::collections::HashMap;

fn offdiag(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1)); for i in 0..n { for j in 0..n { if i != j { v.push(x[i * n + j]); } } } v }
fn upper(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1) / 2); for i in 0..n { for j in i + 1..n { v.push(x[i * n + j]); } } v }
fn sub(full: &[f32], nn: usize, ix: &[usize]) -> Vec<f64> {
    let n = ix.len();
    let mut o = vec![0f64; n * n];
    for (x, &i) in ix.iter().enumerate() { for (y, &j) in ix.iter().enumerate() { if x != y { o[x * n + y] = full[i * nn + j] as f64; } } }
    o
}
/// least squares without intercept for three predictors
fn lsq3(x: [&[f64]; 3], y: &[f64]) -> [f64; 3] {
    let mut a = [[0f64; 4]; 3];
    for p in 0..3 { for q in 0..3 { a[p][q] = x[p].iter().zip(x[q]).map(|(u, v)| u * v).sum(); } a[p][3] = x[p].iter().zip(y).map(|(u, v)| u * v).sum(); }
    for i in 0..3 { for r in 0..3 { if r != i { let f = a[r][i] / a[i][i]; for c in i..4 { a[r][c] -= f * a[i][c]; } } } }
    [a[0][3] / a[0][0], a[1][3] / a[1][1], a[2][3] / a[2][2]]
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 5 { bail!("usage: transfer_prediction <umi_floor> <target> <reference1> <reference2> [...]"); }
    let floor: f64 = a[1].parse()?;
    let paths = &a[2..];
    let scr: Vec<Screen> = paths.iter().map(|p| Screen::load(p)).collect::<Result<_>>()?;
    // genes knocked down (t <= -3) and expressed above the floor in every screen
    let maps: Vec<HashMap<&str, usize>> = scr.iter().map(|d| d.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect()).collect();
    let mut ix: Vec<Vec<usize>> = vec![vec![]; scr.len()];
    for g in &scr[0].genes {
        let pos: Vec<Option<&usize>> = maps.iter().map(|m| m.get(g.as_str())).collect();
        if pos.iter().any(|p| p.is_none()) { continue; }
        let pos: Vec<usize> = pos.into_iter().map(|p| *p.unwrap()).collect();
        if scr.iter().zip(&pos).any(|(d, &i)| d.umi[i] < floor || d.t(&d.mean, i, i, 1.0) > -3.0) { continue; }
        for (k, &i) in pos.iter().enumerate() { ix[k].push(i); }
    }
    let n = ix[0].len();
    let m: Vec<Vec<f64>> = scr.iter().zip(&ix).map(|(d, i)| sub(&d.mean, d.n(), i)).collect();
    let rho: Vec<Vec<f64>> = paths.iter().zip(&scr).zip(&ix).map(|((p, d), i)| Ok(sub(&read_f32(&format!("{p}.ctrlcorr.f32"))?, d.n(), i))).collect::<Result<_>>()?;
    let nref = scr.len() - 1;
    let dot = |x: &[f64], y: &[f64]| -> f64 { x.iter().zip(y).map(|(u, v)| u * v).sum() };
    // scale of co-expression, pooled over the reference lines
    let c = (1..=nref).map(|k| dot(&rho[k], &m[k])).sum::<f64>() / (1..=nref).map(|k| dot(&rho[k], &rho[k])).sum::<f64>();
    // weights fitted between ordered pairs of reference lines
    let mut w = [0f64; 3];
    let mut np = 0f64;
    for x in 1..=nref { for z in 1..=nref { if x != z {
        let v = lsq3([&m[x], &rho[z], &rho[x]], &m[z]);
        for p in 0..3 { w[p] += v[p]; }
        np += 1.0;
    } } }
    for p in 0..3 { w[p] /= np; }

    let mut preds: Vec<(&str, Vec<f64>)> = vec![];
    let mean_ref: Vec<f64> = (0..n * n).map(|i| (1..=nref).map(|k| m[k][i]).sum::<f64>() / nref as f64).collect();
    let mean_rho: Vec<f64> = (0..n * n).map(|i| (1..=nref).map(|k| rho[k][i]).sum::<f64>() / nref as f64).collect();
    let (_, skew_ref) = split(&mean_ref, n);
    preds.push(("target_coexpression", rho[0].iter().map(|r| c * r).collect()));
    preds.push(("reference_response", mean_ref.clone()));
    preds.push(("reference_adjusted", (0..n * n).map(|i| mean_ref[i] + c * (rho[0][i] - mean_rho[i])).collect()));
    preds.push(("reference_skew_target_coexpr", (0..n * n).map(|i| c * rho[0][i] + skew_ref[i]).collect()));
    preds.push(("transfer_fitted", (0..n * n).map(|i| w[0] * mean_ref[i] + w[1] * rho[0][i] + w[2] * mean_rho[i]).collect()));

    let (s0, a0) = split(&m[0], n);
    let (o0, us0, ua0) = (offdiag(&m[0], n), upper(&s0, n), upper(&a0, n));
    println!("target\treferences\tumi_floor\tn_genes\tpredictor\tr_whole\tr_symmetric_part\tr_skew_part\tfitted");
    for (name, p) in &preds {
        let (ps, pa) = split(p, n);
        let ra = if pa.iter().any(|&v| v != 0.0) { pearson(&upper(&pa, n), &ua0) } else { 0.0 };
        let note = match *name { "transfer_fitted" => format!("[{:.3} {:.3} {:.3}]", w[0], w[1], w[2]), "target_coexpression" => format!("{c:.3}"), _ => String::new() };
        println!("{}\t{}\t{floor}\t{n}\t{name}\t{:.4}\t{:.4}\t{ra:.4}\t{note}", paths[0], paths[1..].join(","), pearson(&offdiag(p, n), &o0), pearson(&upper(&ps, n), &us0));
    }
    let (h0, h1) = (sub(&scr[0].h0, scr[0].n(), &ix[0]), sub(&scr[0].h1, scr[0].n(), &ix[0]));
    let ((sa, aa), (sb, ab)) = (split(&h0, n), split(&h1, n));
    println!("{}\t{}\t{floor}\t{n}\ttarget_half_vs_half\t{:.4}\t{:.4}\t{:.4}\t", paths[0], paths[1..].join(","), pearson(&offdiag(&h0, n), &offdiag(&h1, n)), pearson(&upper(&sa, n), &upper(&sb, n)), pearson(&upper(&aa, n), &upper(&ab, n)));
    Ok(())
}
