//! Gene-level view of the skew-symmetric part.
//!
//! Usage: gene_level <screen1> <screen2|SPLIT> <umi_floor> <out_prefix>
//!
//! 1. Additive model: how much of the replicable skew-symmetric variance is a gene-level term
//!    A_ij = p_i - p_j (a knockdown that shifts many genes, a gene shifted by few knockdowns),
//!    and how much is specific to the pair.
//! 2. Sign agreement of the two directions of a pair, which is unaffected by any positive
//!    scaling of rows or columns.
//! 3. For every gene, the replicable effect it sends (as a knockdown) and receives (as a readout).
//! Output: <out_prefix>.summary.tsv, <out_prefix>.genes.tsv
use anyhow::Result;
use perturbseq_reciprocity::{build, pearson, shared_genes, spearman, split, Screen};
use std::fmt::Write as _;

/// residuals of S after s_ij = g + u_i + u_j and of A after a_ij = p_i - p_j (least squares, off-diagonal)
fn residuals(s: &[f64], a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let k = (n - 1) as f64;
    let p: Vec<f64> = (0..n).map(|i| a[i * n..(i + 1) * n].iter().sum::<f64>() / n as f64).collect();
    let rs: Vec<f64> = (0..n).map(|i| s[i * n..(i + 1) * n].iter().sum::<f64>()).collect();
    let g = rs.iter().sum::<f64>() / (n as f64 * k);
    let u: Vec<f64> = rs.iter().map(|r| (r - k * g) / (n as f64 - 2.0)).collect();
    let (mut sr, mut ar) = (s.to_vec(), a.to_vec());
    for i in 0..n { for j in 0..n { if i != j { sr[i * n + j] -= g + u[i] + u[j]; ar[i * n + j] -= p[i] - p[j]; } } }
    (sr, ar, p)
}
fn dot(x: &[f64], y: &[f64]) -> f64 { x.iter().zip(y).map(|(a, b)| a * b).sum() }
fn upper(x: &[f64], n: usize) -> Vec<f64> { let mut v = Vec::with_capacity(n * (n - 1) / 2); for i in 0..n { for j in i + 1..n { v.push(x[i * n + j]); } } v }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let floor: f64 = a[3].parse()?;
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let (ix1, ix2) = shared_genes(&d1, &d2, floor);
    let n = ix1.len();
    let (b1, b2) = (build(&d1, &d1.mean, &ix1, 1.0), build(&d2, &d2.mean, &ix2, 1.0));
    let mut out = String::from("umi_floor\tmatrix\tn_genes\tV_sym\tV_asym\tf_asym\tshare_of_V_asym_gene_level\tf_asym_pair_specific\tr_A_pair_specific\tr_S_pair_specific\tr_gene_level_term\trho_sent_vs_received\n");
    let mut genes = String::from("gene\tmatrix\tgene_level_term_screen1\tgene_level_term_screen2\treplicable_effect_sent\treplicable_effect_received\n");
    for (name, m1, m2) in [("scaled", &b1.sc, &b2.sc), ("raw", &b1.raw, &b2.raw)] {
        let ((s1, a1), (s2, a2)) = (split(m1, n), split(m2, n));
        let ((sr1, ar1, p1), (sr2, ar2, p2)) = (residuals(&s1, &a1, n), residuals(&s2, &a2, n));
        let (vs, va, vsr, var) = (dot(&s1, &s2), dot(&a1, &a2), dot(&sr1, &sr2), dot(&ar1, &ar2));
        // replicable effect sent by knockdown i and received by gene i (cross-screen products of the response)
        let sent: Vec<f64> = (0..n).map(|i| (0..n).filter(|&j| j != i).map(|j| m1[i * n + j] * m2[i * n + j]).sum()).collect();
        let recv: Vec<f64> = (0..n).map(|i| (0..n).filter(|&j| j != i).map(|j| m1[j * n + i] * m2[j * n + i]).sum()).collect();
        writeln!(out, "{floor}\t{name}\t{n}\t{vs:.3}\t{va:.3}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            va / (vs + va), (va - var) / va, var / (vsr + var), pearson(&upper(&ar1, n), &upper(&ar2, n)), pearson(&upper(&sr1, n), &upper(&sr2, n)), pearson(&p1, &p2), spearman(&sent, &recv))?;
        for i in 0..n { writeln!(genes, "{}\t{name}\t{:.5}\t{:.5}\t{:.4}\t{:.4}", d1.genes[ix1[i]], p1[i], p2[i], sent[i], recv[i])?; }
    }
    // sign agreement of the two directions; pairs chosen in one screen and read in the other
    out.push_str("select\tread\tn_pairs_both_directions_t3\tsame_sign\tfraction\n");
    for (sel, ev, ls, le) in [(&b1, &b2, "screen1", "screen2"), (&b2, &b1, "screen2", "screen1")] {
        let (mut c, mut k) = (0usize, 0usize);
        for i in 0..n { for j in i + 1..n {
            if sel.t[i * n + j].abs() >= 3.0 && sel.t[j * n + i].abs() >= 3.0 { c += 1; if ev.raw[i * n + j] * ev.raw[j * n + i] > 0.0 { k += 1; } }
        } }
        writeln!(out, "{ls}\t{le}\t{c}\t{k}\t{:.4}", k as f64 / c.max(1) as f64)?;
    }
    std::fs::write(format!("{}.summary.tsv", a[4]), out)?;
    std::fs::write(format!("{}.genes.tsv", a[4]), genes)?;
    Ok(())
}
