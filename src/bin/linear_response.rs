//! Linear-response prediction from co-expression in control cells, and its ceiling.
//!
//! Usage: linear_response <screen1> <screen2|SPLIT> <umi_floor>   (TSV on standard output)
//!
//! On the standardised scale the prediction for knockdown A is Mhat_AB = u_A * rho_AB, where
//! rho is the correlation of standardised expression in the control cells of screen 1.
//!   fitted scale    : u_A is the single free parameter, fitted by least squares on screen 1
//!   on-target scale : u_A is the on-target shift of A, no free parameter
//!   symmetric       : one global scale, Mhat = c * rho
//! Each prediction is compared with the independent screen 2. The agreement of screen 1 with
//! screen 2 is the ceiling set by measurement noise.
use anyhow::Result;
use perturbseq_reciprocity::{build, pearson, read_f32, shared_genes, spearman, split, Screen};

fn offdiag(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1)); for i in 0..n { for j in 0..n { if i != j { v.push(x[i * n + j]); } } } v }
fn upper(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1) / 2); for i in 0..n { for j in i + 1..n { v.push(x[i * n + j]); } } v }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let floor: f64 = a[3].parse()?;
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let rho_full = read_f32(&format!("{}.ctrlcorr.f32", a[1]))?;
    let (ix1, ix2) = shared_genes(&d1, &d2, floor);
    let n = ix1.len();
    let nn = d1.n();
    let (m1, m2) = (build(&d1, &d1.mean, &ix1, 1.0).raw, build(&d2, &d2.mean, &ix2, 1.0).raw);
    let mut rho = vec![0f64; n * n];
    for (x, &ia) in ix1.iter().enumerate() { for (y, &ib) in ix1.iter().enumerate() { rho[x * n + y] = rho_full[ia * nn + ib] as f64; } }
    let (mut fit, mut strict, mut sym) = (vec![0f64; n * n], vec![0f64; n * n], vec![0f64; n * n]);
    let (mut gnum, mut gden) = (0f64, 0f64);
    for i in 0..n { for j in 0..n { if i != j { gnum += rho[i * n + j] * m1[i * n + j]; gden += rho[i * n + j] * rho[i * n + j]; } } }
    for i in 0..n {
        let (mut num, mut den) = (0f64, 0f64);
        for j in 0..n { if i != j { num += rho[i * n + j] * m1[i * n + j]; den += rho[i * n + j] * rho[i * n + j]; } }
        for j in 0..n { if i != j {
            fit[i * n + j] = num / den * rho[i * n + j];
            strict[i * n + j] = m1[i * n + i] * rho[i * n + j];
            sym[i * n + j] = gnum / gden * rho[i * n + j];
        } }
    }
    let ((s1, a1), (s2, a2)) = (split(&m1, n), split(&m2, n));
    let (o2, us2, ua2) = (offdiag(&m2, n), upper(&s2, n), upper(&a2, n));
    let ceil = pearson(&offdiag(&m1, n), &o2);
    let (cs, ca) = (pearson(&upper(&s1, n), &us2), pearson(&upper(&a1, n), &ua2));
    let (vs, va): (f64, f64) = (s1.iter().zip(&s2).map(|(x, y)| x * y).sum(), a1.iter().zip(&a2).map(|(x, y)| x * y).sum());
    println!("screen1\tscreen2\tumi_floor\tpredictor\tn_genes\tr_prediction\tr_ceiling\tr_symmetric_part\tr_ceiling_symmetric\tr_skew_part\tr_ceiling_skew\tskew_share_of_prediction\tf_asym_observed\trho_error_vs_abs_skew\trho_error_vs_abs_symmetric");
    for (name, pred) in [("fitted_scale", &fit), ("on_target_scale", &strict), ("symmetric", &sym)] {
        let (ps, pa) = split(pred, n);
        let op = offdiag(pred, n);
        let (pvs, pva): (f64, f64) = (ps.iter().map(|x| x * x).sum(), pa.iter().map(|x| x * x).sum());
        let slope = op.iter().zip(&o2).map(|(x, y)| x * y).sum::<f64>() / op.iter().map(|x| x * x).sum::<f64>();
        let err: Vec<f64> = op.iter().zip(&o2).map(|(x, y)| (y - slope * x).powi(2)).collect();
        let (aa, ss): (Vec<f64>, Vec<f64>) = (offdiag(&a1, n).iter().map(|v| v.abs()).collect(), offdiag(&s1, n).iter().map(|v| v.abs()).collect());
        let ra = if pva > 0.0 { pearson(&upper(&pa, n), &ua2) } else { 0.0 };
        println!("{}\t{}\t{floor}\t{name}\t{n}\t{:.4}\t{ceil:.4}\t{:.4}\t{cs:.4}\t{ra:.4}\t{ca:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            a[1], a[2], pearson(&op, &o2), pearson(&upper(&ps, n), &us2), pva / (pvs + pva), va / (vs + va), spearman(&err, &aa), spearman(&err, &ss));
    }
    Ok(())
}
