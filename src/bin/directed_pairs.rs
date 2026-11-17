//! Counts of gene pairs with a replicated effect in one direction or in both.
//!
//! Usage: directed_pairs <screen1> <screen2|SPLIT> <umi_floor> <out_prefix>
//!
//! An effect of knockdown i on gene j is called replicated when |t| >= 3 in both screens with
//! the same sign. For every pair with at least one replicated direction the reverse direction
//! is classified as replicated, absent (|t| < 2 in both screens) or intermediate.
//! Output: <out_prefix>.counts.tsv and <out_prefix>.pairs.tsv (replicated one-direction pairs
//! with an absent reverse direction).
use anyhow::Result;
use perturbseq_reciprocity::{build, shared_genes, Screen};
use std::fmt::Write as _;

fn median(v: &mut Vec<f64>) -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v[v.len() / 2] }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let floor: f64 = a[3].parse()?;
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let (ix1, ix2) = shared_genes(&d1, &d2, floor);
    let n = ix1.len();
    let (b1, b2) = (build(&d1, &d1.mean, &ix1, 1.0), build(&d2, &d2.mean, &ix2, 1.0));
    let rep = |i: usize, j: usize| { let (x, y) = (b1.t[i * n + j], b2.t[i * n + j]); x.abs() >= 3.0 && y.abs() >= 3.0 && x * y > 0.0 };
    let absent = |i: usize, j: usize| b1.t[i * n + j].abs() < 2.0 && b2.t[i * n + j].abs() < 2.0;
    let eff = |i: usize, j: usize| 0.5 * (b1.sc[i * n + j] + b2.sc[i * n + j]);
    let (mut both, mut same, mut one_absent, mut one_mid, mut none) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut ratio, mut ratio_both) = (vec![], vec![]);
    let mut pairs = String::from("knocked_down_gene\tresponding_gene\tt_screen1\tt_screen2\tscaled_response\treverse_t_screen1\treverse_t_screen2\treverse_scaled_response\n");
    for i in 0..n {
        for j in i + 1..n {
            let (f, r) = (rep(i, j), rep(j, i));
            if f && r {
                both += 1;
                if eff(i, j) * eff(j, i) > 0.0 { same += 1; }
                let (x, y) = (eff(i, j).abs(), eff(j, i).abs());
                ratio_both.push(x.min(y) / x.max(y));
            } else if f || r {
                let (u, v) = if f { (i, j) } else { (j, i) };
                ratio.push(eff(v, u).abs() / eff(u, v).abs());
                if absent(v, u) {
                    one_absent += 1;
                    writeln!(pairs, "{}\t{}\t{:.2}\t{:.2}\t{:.4}\t{:.2}\t{:.2}\t{:.4}", d1.genes[ix1[u]], d1.genes[ix1[v]], b1.t[u * n + v], b2.t[u * n + v], eff(u, v), b1.t[v * n + u], b2.t[v * n + u], eff(v, u))?;
                } else { one_mid += 1; }
            } else { none += 1; }
        }
    }
    let with = both + one_absent + one_mid;
    let mut out = String::from("screen1\tscreen2\tumi_floor\tn_genes\tn_pairs\tpairs_with_replicated_effect\tboth_directions\tboth_directions_same_sign\tone_direction_reverse_absent\tone_direction_reverse_intermediate\tfraction_both_directions\tmedian_reverse_to_forward_ratio_one_direction\tmedian_smaller_to_larger_ratio_both_directions\n");
    writeln!(out, "{}\t{}\t{floor}\t{n}\t{}\t{with}\t{both}\t{same}\t{one_absent}\t{one_mid}\t{:.4}\t{:.4}\t{:.4}",
        a[1], a[2], with + none, both as f64 / with.max(1) as f64, median(&mut ratio), median(&mut ratio_both))?;
    std::fs::write(format!("{}.counts.tsv", a[4]), out)?;
    std::fs::write(format!("{}.pairs.tsv", a[4]), pairs)?;
    Ok(())
}
