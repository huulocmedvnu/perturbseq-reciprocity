//! Prediction of the response to knockdowns that are held out of the training screen.
//!
//! Usage: heldout_prediction <screen1> <screen2|SPLIT> <umi_floor> <tf_target.tsv> <complexes.tsv> <folds> [raw1 [raw2]]
//!        (TSV on standard output)
//! With raw1 (and raw2 for a second screen), written by `reduce <cache> <out> raw`, the response that is
//! predicted and scored is the shift of the mean unnormalised count, and the co-expression baselines use
//! the covariance of unnormalised counts in control cells, the formulation of Kuznets-Speck et al.
//! The gene set, the folds and the inputs of the ridge models stay those of the standardised screens.
//!
//! The knocked-down genes are divided into folds. For a held-out gene A the whole row of the
//! response matrix of screen 1 (A knocked down, every gene measured) is hidden. What remains
//! known about A is its correlation with every gene in control cells, its expression after
//! the knockdown of every training gene (the column of A restricted to training rows), and
//! curated annotations. Every predictor is fitted on screen 1 and compared with the rows of
//! the independent screen 2, for the whole response and for its symmetric and skew-symmetric
//! parts. The agreement of screen 1 with screen 2 is the ceiling set by measurement noise.
//!
//!   coexpression_global_scale     c * rho_AB, one scale c fitted on the training knockdowns
//!   coexpression_on_target_scale  (on-target shift of A) * rho_AB
//!   training_mean                 mean response of gene B over the training knockdowns
//!   reverse_direction             response of A to the knockdown of B, for B in training
//!   ridge_coexpression            kernel ridge on the co-expression profile of A
//!   ridge_received                kernel ridge on the received-response profile of A
//!   ridge_received_coexpression   both profiles
//!   ridge_annotations             kernel ridge on annotations of A (regulators, targets, complexes)
//!   ridge_all                     both profiles and annotations
//!   stacked_simple                least-squares combination of training mean, co-expression and reverse direction
//!   stacked_all                   the same with the three single-source ridge models added
//!   response_scale_from_screen1   u_A * rho_AB, u_A fitted to the response to knockdown A in screen 1
//!   response_scale_from_screen2   u_A * rho_AB, u_A fitted to the response in screen 2, the response that is scored
//!   mean_expression_global_scale  c * (mean count of B in control cells), one scale fitted on the training knockdowns
//!   mean_expression_scale_from_screen2  the same with one scale per knockdown fitted to the response in screen 2
//!   reference_lines               mean response to the knockdown of A in screens of other cell lines (refs=...)
//!   stacked_with_reference_lines  stacked_all with the reference response added
//! The predictors with a scale from screen 1 or screen 2 use the response to the knockdown that is predicted and are not held-out predictions.
//! The ridge models predict the row of A after subtraction of the training mean. The penalty is
//! chosen by exact leave-one-out prediction within the training knockdowns. stacked_all combines
//! training mean, co-expression, reverse direction and the three single-source ridge models with
//! least-squares weights fitted to their leave-one-out predictions on the training knockdowns.
use anyhow::Result;
use perturbseq_reciprocity::{build, pearson, read_f32, shared_genes, split, Screen};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;

const LAMBDA: [f64; 8] = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0];

fn offdiag(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1)); for i in 0..n { for j in 0..n { if i != j { v.push(x[i * n + j]); } } } v }
fn upper(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1) / 2); for i in 0..n { for j in i + 1..n { v.push(x[i * n + j]); } } v }

/// Gram matrix of the rows of x (n rows, d columns), divided by its mean diagonal
fn gram(x: &[f64], n: usize, d: usize) -> Vec<f64> {
    let mut k = vec![0f64; n * n];
    k.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        let xi = &x[i * d..(i + 1) * d];
        for j in 0..n { row[j] = xi.iter().zip(&x[j * d..(j + 1) * d]).map(|(a, b)| a * b).sum(); }
    });
    let tr = (0..n).map(|i| k[i * n + i]).sum::<f64>() / n as f64;
    if tr > 0.0 { for v in k.iter_mut() { *v /= tr; } }
    k
}

/// Kernel ridge through the eigendecomposition of the training kernel. Rows `tr` are fitted and
/// rows `te` predicted, y holds one row of m targets per training row. For every penalty the exact
/// leave-one-out predictions of the training rows follow in closed form. The penalty with the
/// highest leave-one-out correlation is used. Returns (penalty, leave-one-out predictions, test predictions).
fn ridge(k: &[f64], n: usize, tr: &[usize], te: &[usize], y: &[f64], m: usize) -> (f64, Vec<f64>, Vec<f64>) {
    let t = tr.len();
    let e = nalgebra::DMatrix::from_fn(t, t, |i, j| k[tr[i] * n + tr[j]]).symmetric_eigen();
    let w: Vec<f64> = e.eigenvalues.iter().map(|&x| x.max(0.0)).collect();
    // q[i * t + c] = component i of eigenvector c; qt is its transpose
    let (mut q, mut qt) = (vec![0f64; t * t], vec![0f64; t * t]);
    for i in 0..t { for c in 0..t { let v = e.eigenvectors[(i, c)]; q[i * t + c] = v; qt[c * t + i] = v; } }
    let mul = |a: &[f64], b: &[f64], scale: Option<&[f64]>| -> Vec<f64> { // (t x t) times (t x m), rows of b optionally scaled
        let mut o = vec![0f64; t * m];
        o.par_chunks_mut(m).enumerate().for_each(|(i, row)| {
            for c in 0..t {
                let f = a[i * t + c] * scale.map_or(1.0, |s| s[c]);
                if f != 0.0 { for (x, v) in row.iter_mut().zip(&b[c * m..(c + 1) * m]) { *x += f * v; } }
            }
        });
        o
    };
    let z = mul(&qt, y, None);
    let mut best = (f64::MIN, 0.0, Vec::new(), Vec::new());
    for &lam in &LAMBDA {
        let d: Vec<f64> = w.iter().map(|x| 1.0 / (x + lam)).collect();
        let alpha = mul(&q, &z, Some(&d));
        let hd: Vec<f64> = (0..t).map(|i| (0..t).map(|c| q[i * t + c] * q[i * t + c] * d[c]).sum()).collect();
        let mut loo = vec![0f64; t * m];
        loo.par_chunks_mut(m).enumerate().for_each(|(i, row)| { for b in 0..m { row[b] = y[i * m + b] - alpha[i * m + b] / hd[i]; } });
        let (mut u, mut v) = (Vec::with_capacity(t * m), Vec::with_capacity(t * m));
        for (i, &g) in tr.iter().enumerate() { for b in 0..m { if b != g { u.push(loo[i * m + b]); v.push(y[i * m + b]); } } }
        let r = pearson(&u, &v);
        if r > best.0 { best = (r, lam, loo, alpha); }
    }
    let alpha = &best.3;
    let mut out = vec![0f64; te.len() * m];
    out.par_chunks_mut(m).zip(te.par_iter()).for_each(|(row, &g)| {
        for (i, &h) in tr.iter().enumerate() {
            let f = k[g * n + h];
            if f != 0.0 { for (o, v) in row.iter_mut().zip(&alpha[i * m..(i + 1) * m]) { *o += f * v; } }
        }
    });
    (best.1, best.2, out)
}

/// least squares without intercept, normal equations solved by Gaussian elimination
fn lsq(xtx: &[f64], xty: &[f64], p: usize) -> Vec<f64> {
    let mut a: Vec<Vec<f64>> = (0..p).map(|i| { let mut r = xtx[i * p..(i + 1) * p].to_vec(); r.push(xty[i]); r }).collect();
    for i in 0..p {
        let piv = (i..p).max_by(|&x, &y| a[x][i].abs().partial_cmp(&a[y][i].abs()).unwrap()).unwrap();
        a.swap(i, piv);
        if a[i][i].abs() < 1e-12 { continue; }
        for r in 0..p { if r != i { let f = a[r][i] / a[i][i]; for c in i..=p { a[r][c] -= f * a[i][c]; } } }
    }
    (0..p).map(|i| if a[i][i].abs() < 1e-12 { 0.0 } else { a[i][p] / a[i][i] }).collect()
}

fn main() -> Result<()> {
    let all: Vec<String> = std::env::args().collect();
    // refs=<screen>,<screen>: screens of other cell lines in which the held-out knockdowns were measured
    let refs: Vec<String> = all.iter().find_map(|x| x.strip_prefix("refs=")).map(|x| x.split(',').map(String::from).collect()).unwrap_or_default();
    let a: Vec<String> = all.into_iter().filter(|x| !x.starts_with("refs=")).collect();
    let floor: f64 = a[3].parse()?;
    let nfold: usize = a.get(6).map(|s| s.parse().unwrap()).unwrap_or(5);
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let rho_full = read_f32(&format!("{}.ctrlcorr.f32", a[1]))?;
    let (ix1, ix2) = shared_genes(&d1, &d2, floor);
    let n = ix1.len();
    let nn = d1.n();
    let m1 = build(&d1, &d1.mean, &ix1, 1.0).raw;
    let mut rho = vec![0f64; n * n];
    for (x, &ia) in ix1.iter().enumerate() { for (y, &ib) in ix1.iter().enumerate() { if x != y { rho[x * n + y] = rho_full[ia * nn + ib] as f64; } } }
    // y1, y2: the response that is fitted (screen 1) and scored (screen 2); sig: co-expression of the baselines
    let (y1, y2, sig, sig_diag, scale) = if let Some(r1) = a.get(7) {
        let (e1, e2) = Screen::pair(r1, if a[2] == "SPLIT" { "SPLIT" } else { a.get(8).expect("raw2") })?;
        let pick = |e: &Screen, d: &Screen, ix: &[usize]| -> Vec<f64> {
            let p: HashMap<&str, usize> = e.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
            let j: Vec<usize> = ix.iter().map(|&i| p[d.genes[i].as_str()]).collect();
            let mut o = vec![0f64; n * n];
            for x in 0..n { for y in 0..n { o[x * n + y] = e.mean[j[x] * e.n() + j[y]] as f64; } }
            o
        };
        let (cr, sd) = (read_f32(&format!("{r1}.ctrlcorr.f32"))?, read_f32(&format!("{r1}.ctrlsd.f32"))?);
        let p: HashMap<&str, usize> = e1.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
        let j: Vec<usize> = ix1.iter().map(|&i| p[d1.genes[i].as_str()]).collect();
        let mut sg = vec![0f64; n * n];
        for x in 0..n { for y in 0..n { if x != y { sg[x * n + y] = cr[j[x] * e1.n() + j[y]] as f64 * sd[j[x]] as f64 * sd[j[y]] as f64; } } }
        let dg: Vec<f64> = j.iter().map(|&i| (sd[i] as f64).powi(2)).collect();
        (pick(&e1, &d1, &ix1), pick(&e2, &d2, &ix2), sg, dg, "unnormalised_counts")
    } else {
        (m1.clone(), build(&d2, &d2.mean, &ix2, 1.0).raw, rho.clone(), vec![1.0; n], "standardised")
    };
    // mean response to the same knockdown in the reference lines (standardised scale), where measured
    let mut refm = vec![0f64; n * n];
    let mut refc = vec![0f64; n * n];
    for r in &refs {
        let e = Screen::load(r)?;
        let p: HashMap<&str, usize> = e.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
        let j: Vec<Option<usize>> = ix1.iter().map(|&i| p.get(d1.genes[i].as_str()).copied().filter(|&q| e.t(&e.mean, q, q, 1.0) <= -3.0)).collect();
        for x in 0..n { let Some(jx) = j[x] else { continue }; for y in 0..n { if let Some(jy) = p.get(d1.genes[ix1[y]].as_str()) { refm[x * n + y] += e.mean[jx * e.n() + jy] as f64; refc[x * n + y] += 1.0; } } }
    }
    for i in 0..n * n { if refc[i] > 0.0 { refm[i] /= refc[i]; } }
    let expr: Vec<f64> = ix1.iter().map(|&i| d1.umi[i]).collect();
    let pos: HashMap<&str, usize> = ix1.iter().enumerate().map(|(k, &i)| (d1.genes[i].as_str(), k)).collect();

    // annotations of gene k: its targets, its regulators, the members of its complexes
    let mut ann = vec![0f64; n * 3 * n];
    for l in std::fs::read_to_string(&a[4])?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 2 { continue; }
        if let (Some(&t), Some(&g)) = (pos.get(f[0]), pos.get(f[1])) { if t != g { ann[t * 3 * n + g] = 1.0; ann[g * 3 * n + n + t] = 1.0; } }
    }
    for l in std::fs::read_to_string(&a[5])?.lines() {
        let Some(m) = l.split('\t').nth(1) else { continue };
        let g: Vec<usize> = m.split('_').filter_map(|x| pos.get(x).copied()).collect();
        for &x in &g { for &y in &g { if x != y { ann[x * 3 * n + 2 * n + y] = 1.0; } } }
    }
    let k_rho = gram(&rho, n, n);
    let k_ann = gram(&ann, n, 3 * n);
    drop(ann);

    // folds
    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(&mut Xoshiro256PlusPlus::seed_from_u64(11));
    let mut fold = vec![0usize; n];
    for (p, &g) in order.iter().enumerate() { fold[g] = p % nfold; }

    let names = ["coexpression_global_scale", "coexpression_on_target_scale", "training_mean", "reverse_direction", "ridge_coexpression",
        "ridge_received", "ridge_received_coexpression", "ridge_all", "stacked_simple", "stacked_all", "ridge_annotations",
        "response_scale_from_screen1", "response_scale_from_screen2", "mean_expression_global_scale", "mean_expression_scale_from_screen2",
        "reference_lines", "stacked_with_reference_lines"];
    let mut pred: Vec<Vec<f64>> = names.iter().map(|_| vec![0f64; n * n]).collect();
    let mut notes: Vec<String> = vec![String::new(); names.len()];

    for f in 0..nfold {
        let tr: Vec<usize> = (0..n).filter(|&g| fold[g] != f).collect();
        let te: Vec<usize> = (0..n).filter(|&g| fold[g] == f).collect();
        let is_tr: Vec<bool> = (0..n).map(|g| fold[g] != f).collect();
        // one global scale of co-expression, and the mean response of every gene, over training rows
        let (mut num, mut den) = (0f64, 0f64);
        let mut mean = vec![0f64; n];
        let mut cnt = vec![0f64; n];
        for &k in &tr { for b in 0..n { if b != k { num += sig[k * n + b] * y1[k * n + b]; den += sig[k * n + b] * sig[k * n + b]; mean[b] += y1[k * n + b]; cnt[b] += 1.0; } } }
        let c = num / den;
        let (mut nu, mut de) = (0f64, 0f64);
        for &k in &tr { for b in 0..n { if b != k { nu += expr[b] * y1[k * n + b]; de += expr[b] * expr[b]; } } }
        let ce = nu / de;
        for b in 0..n { mean[b] /= cnt[b]; }
        // received-response profile of every gene: its shift under each training knockdown
        let t = tr.len();
        let mut recv = vec![0f64; n * t];
        for g in 0..n { for (j, &k) in tr.iter().enumerate() { if k != g { recv[g * t + j] = m1[k * n + g]; } } }
        let k_recv = gram(&recv, n, t);
        drop(recv);
        let add = |x: &[f64], y: &[f64]| -> Vec<f64> { x.iter().zip(y).map(|(u, v)| u + v).collect() };
        let k_rc = add(&k_recv, &k_rho);
        let k_all = add(&k_rc, &k_ann);
        // targets: training rows minus the training mean, on-target entry set to zero
        let mut y = vec![0f64; t * n];
        for (i, &k) in tr.iter().enumerate() { for b in 0..n { if b != k { y[i * n + b] = y1[k * n + b] - mean[b]; } } }
        let kernels: [(&[f64], usize); 5] = [(&k_rho, 4), (&k_recv, 5), (&k_rc, 6), (&k_all, 7), (&k_ann, 10)];
        // leave-one-out predictions of the three single-source models on the training rows, and their test predictions
        let mut loo: Vec<Vec<f64>> = vec![vec![0f64; t * n]; 3];
        let mut tst: Vec<Vec<f64>> = vec![vec![0f64; te.len() * n]; 3];
        for (kern, slot) in kernels {
            if std::env::var("BASELINES_ONLY").is_ok() { continue; } // ridge models skipped, their rows are then not meaningful
            let (lam, l, p) = ridge(kern, n, &tr, &te, &y, n);
            notes[slot] += &format!("{}{lam}", if f > 0 { "," } else { "" });
            for (r, &g) in te.iter().enumerate() { for b in 0..n { if b != g { pred[slot][g * n + b] = mean[b] + p[r * n + b]; } } }
            if let Some(x) = [4usize, 5, 10].iter().position(|&v| v == slot) { loo[x] = l; tst[x] = p; }
        }
        // stacking weights on the training rows; the reverse direction is known only for training genes B
        let base = |k: usize, b: usize| -> [f64; 3] { [mean[b], c * sig[k * n + b], if is_tr[b] { y1[b * n + k] } else { 0.0 }] };
        let (mut xtx, mut xty) = (vec![0f64; 49], vec![0f64; 7]);
        for (i, &k) in tr.iter().enumerate() { for b in 0..n { if b != k {
            let h = base(k, b);
            let x = [h[0], h[1], h[2], loo[0][i * n + b], loo[1][i * n + b], loo[2][i * n + b], refm[k * n + b]];
            for p in 0..7 { xty[p] += x[p] * y1[k * n + b]; for q in 0..7 { xtx[p * 7 + q] += x[p] * x[q]; } }
        } } }
        let sub = |m: usize| -> Vec<f64> { lsq(&(0..m).flat_map(|p| (0..m).map(move |q| (p, q))).map(|(p, q)| xtx[p * 7 + q]).collect::<Vec<_>>(), &xty[..m], m) };
        let (w3, w6) = (sub(3), sub(6));
        let w7 = if refs.is_empty() { vec![0.0; 7] } else { sub(7) };
        // the reference response alone, on the scale of this screen
        let wr = xty[6] / xtx[48].max(1e-300);
        let show = |w: &[f64]| -> String { format!("{}[{}]", if f > 0 { "," } else { "" }, w.iter().map(|x| format!("{x:.3}")).collect::<Vec<_>>().join(" ")) };
        notes[8] += &show(&w3);
        notes[9] += &show(&w6);
        if !refs.is_empty() { notes[16] += &show(&w7); }
        for (r, &g) in te.iter().enumerate() { for b in 0..n { if b != g {
            let i = g * n + b;
            pred[0][i] = c * sig[i];
            pred[1][i] = y1[g * n + g] / sig_diag[g] * sig[i];
            pred[2][i] = mean[b];
            pred[13][i] = ce * expr[b];
            pred[3][i] = if is_tr[b] { y1[b * n + g] } else { c * sig[i] };
            let h = base(g, b);
            let x = [h[0], h[1], h[2], tst[0][r * n + b], tst[1][r * n + b], tst[2][r * n + b], refm[i]];
            pred[8][i] = (0..3).map(|p| w3[p] * x[p]).sum();
            pred[9][i] = (0..6).map(|p| w6[p] * x[p]).sum();
            pred[15][i] = wr * refm[i];
            pred[16][i] = (0..7).map(|p| w7[p] * x[p]).sum();
        } } }
        eprintln!("fold {} of {nfold}: {} training and {} held-out knockdowns", f + 1, tr.len(), te.len());
    }

    // scale of each knockdown fitted to its own response: in screen 1, and in the scored screen 2
    for g in 0..n {
        let row = |y: &[f64]| -> f64 { let (mut u, mut v) = (0f64, 0f64); for b in 0..n { if b != g { u += sig[g * n + b] * y[g * n + b]; v += sig[g * n + b] * sig[g * n + b]; } } u / v };
        let (u1, u2) = (row(&y1), row(&y2));
        let (mut nu, mut de) = (0f64, 0f64);
        for b in 0..n { if b != g { nu += expr[b] * y2[g * n + b]; de += expr[b] * expr[b]; } }
        for b in 0..n { if b != g { pred[11][g * n + b] = u1 * sig[g * n + b]; pred[12][g * n + b] = u2 * sig[g * n + b]; pred[14][g * n + b] = nu / de * expr[b]; } }
    }
    let (s2, a2) = split(&y2, n);
    let (o2, us2, ua2) = (offdiag(&y2, n), upper(&s2, n), upper(&a2, n));
    let median = |mut v: Vec<f64>| -> f64 { v.retain(|x| x.is_finite()); v.sort_by(|a, b| a.partial_cmp(b).unwrap()); if v.is_empty() { f64::NAN } else { v[v.len() / 2] } };
    // agreement over the genes of one knockdown, the knocked-down gene left out
    let per_row = |p: &[f64]| -> Vec<f64> { (0..n).into_par_iter().map(|g| {
        let (u, v): (Vec<f64>, Vec<f64>) = (0..n).filter(|&b| b != g).map(|b| (p[g * n + b], y2[g * n + b])).unzip();
        pearson(&u, &v) }).collect() };
    println!("screen1\tscreen2\tresponse\tumi_floor\tn_genes\tfolds\tpredictor\tr_whole\tr_symmetric_part\tr_skew_part\tskew_share_of_prediction\tmedian_r_per_knockdown\tmedian_r2_per_knockdown\tmedian_norm_score_per_knockdown\tfitted");
    let line = |name: &str, p: &[f64], note: &str| {
        let (ps, pa) = split(p, n);
        let (vs, va): (f64, f64) = (ps.iter().map(|x| x * x).sum(), pa.iter().map(|x| x * x).sum());
        let ra = if va > 0.0 { pearson(&upper(&pa, n), &ua2) } else { 0.0 };
        let rr = per_row(p);
        // 1 - |measured - predicted| / |measured| over the genes of one knockdown (Kuznets-Speck et al., their R2)
        let ns: Vec<f64> = (0..n).into_par_iter().map(|g| {
            let (mut e, mut t) = (0f64, 0f64);
            for b in 0..n { if b != g { e += (y2[g * n + b] - p[g * n + b]).powi(2); t += y2[g * n + b].powi(2); } }
            1.0 - (e / t).sqrt() }).collect();
        println!("{}\t{}\t{scale}\t{floor}\t{n}\t{nfold}\t{name}\t{:.4}\t{:.4}\t{ra:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{note}", a[1], a[2], pearson(&offdiag(p, n), &o2), pearson(&upper(&ps, n), &us2), va / (vs + va),
            median(rr.clone()), median(rr.iter().map(|r| r * r).collect()), median(ns));
    };
    for (i, name) in names.iter().enumerate() { if i < 15 || !refs.is_empty() { line(name, &pred[i], &notes[i]); } }
    line("ceiling_screen1", &y1, "");
    // the co-expression row with the knocked-down gene included and the scale fitted to the scored response,
    // as in Kuznets-Speck et al.: squared correlation, and 1 - |residual| / |measured|
    let with_target: Vec<(f64, f64)> = (0..n).into_par_iter().map(|g| {
        let (u, v): (Vec<f64>, Vec<f64>) = (0..n).map(|b| (if b == g { sig_diag[g] } else { sig[g * n + b] }, y2[g * n + b])).unzip();
        let (uv, uu, vv): (f64, f64, f64) = (u.iter().zip(&v).map(|(x, y)| x * y).sum(), u.iter().map(|x| x * x).sum(), v.iter().map(|x| x * x).sum());
        (pearson(&u, &v).powi(2), 1.0 - (1.0 - uv * uv / (uu * vv)).max(0.0).sqrt()) }).collect();
    println!("{}\t{}\t{scale}\t{floor}\t{n}\t{nfold}\tcoexpression_scale_from_screen2_with_target_gene\t\t\t\t\t\t{:.4}\t{:.4}\t", a[1], a[2],
        median(with_target.iter().map(|x| x.0).collect()), median(with_target.iter().map(|x| x.1).collect()));
    Ok(())
}
