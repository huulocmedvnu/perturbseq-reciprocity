//! Reciprocity of the gene-by-gene response matrix across two independent screens.
//!
//! Usage: reciprocity <screen1> <screen2|SPLIT> <tf_target.tsv> <complexes.tsv> <tss.tsv> <out_prefix> [n_perm] [n_boot]
//!
//! screen1 supplies the discovery data (active pairs are selected in one half of its cells
//! and read in the other half), screen2 is the replication screen. With SPLIT the two
//! screens are disjoint cell quarters of screen1.
//! Output: <out_prefix>.json (all estimates), <out_prefix>.density.tsv (binned scatter of the
//! symmetric and skew-symmetric parts in the two screens, for plotting).
use anyhow::{bail, Result};
use perturbseq_reciprocity::{binom_p, build, shared_genes, spearman, split, Mat, Screen};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::time::Instant;

#[derive(Clone, Copy, Default)]
struct Acc { n: f64, x: f64, y: f64, xx: f64, yy: f64, xy: f64 }
impl Acc {
    #[inline]
    fn add(&mut self, x: f64, y: f64) { self.n += 1.0; self.x += x; self.y += y; self.xx += x * x; self.yy += y * y; self.xy += x * y; }
    fn r(&self) -> f64 {
        let c = self.xy - self.x * self.y / self.n;
        let vx = self.xx - self.x * self.x / self.n;
        let vy = self.yy - self.y * self.y / self.n;
        c / (vx * vy).sqrt()
    }
}
#[derive(Clone, Copy)]
struct Stat { npair: f64, r_s: f64, r_a: f64, v_sym: f64, v_asym: f64, f_asym: f64 }
struct Pairs<'a> { n: usize, s1: &'a [f64], a1: &'a [f64], s2: &'a [f64], a2: &'a [f64], mask: &'a [u8], need: u8, forbid: u8 }
impl Pairs<'_> {
    #[inline]
    fn keep(&self, u: usize, v: usize) -> bool { let m = self.mask[u * self.n + v]; m & self.need == self.need && m & self.forbid == 0 }
    /// i1 indexes screen 1 (and the pair classes), i2 indexes screen 2
    fn stat(&self, i1: &[usize], i2: &[usize]) -> Stat {
        let n = self.n;
        let (mut s, mut a) = (Acc::default(), Acc::default());
        for p in 0..i1.len() {
            let (u, u2) = (i1[p], i2[p]);
            for q in p + 1..i1.len() {
                let (v, v2) = (i1[q], i2[q]);
                if u == v || !self.keep(u, v) { continue; }
                s.add(self.s1[u * n + v], self.s2[u2 * n + v2]);
                a.add(self.a1[u * n + v], self.a2[u2 * n + v2]);
            }
        }
        let (vs, va) = (2.0 * s.xy, 2.0 * a.xy);
        Stat { npair: s.n, r_s: s.r(), r_a: a.r(), v_sym: vs, v_asym: va, f_asym: va / (vs + va) }
    }
}
fn quant(v: &mut Vec<f64>, q: f64) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * q).round() as usize]
}

const ACTIVE: u8 = 1; // at least one direction with |t| >= 3 in the discovery half of screen 1
const TF: u8 = 2; // curated transcription factor and target
const CPX: u8 = 4; // subunits of one protein complex
const CIS: u8 = 8; // transcription start sites on one chromosome within 1 Mb

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 7 { bail!("usage: reciprocity <screen1> <screen2|SPLIT> <tf_target.tsv> <complexes.tsv> <tss.tsv> <out_prefix> [n_perm] [n_boot]"); }
    let nperm: usize = a.get(7).map(|s| s.parse().unwrap()).unwrap_or(10_000);
    let nboot: usize = a.get(8).map(|s| s.parse().unwrap()).unwrap_or(2_000);
    let t0 = Instant::now();
    let (d1, d2) = Screen::pair(&a[1], &a[2])?;
    let tf: HashSet<(String, String)> = std::fs::read_to_string(&a[3])?.lines()
        .filter_map(|l| { let f: Vec<&str> = l.split('\t').collect(); (f.len() >= 2 && f[0] != f[1]).then(|| (f[0].to_string(), f[1].to_string())) }).collect();
    let complexes: Vec<Vec<String>> = std::fs::read_to_string(&a[4])?.lines()
        .filter_map(|l| l.split('\t').nth(1).map(|m| m.split('_').map(String::from).collect())).collect();
    let tss: HashMap<String, (String, f64)> = std::fs::read_to_string(&a[5])?.lines()
        .filter_map(|l| { let f: Vec<&str> = l.split('\t').collect(); (f.len() >= 3).then(|| (f[0].to_string(), (f[1].to_string(), f[2].parse().unwrap_or(0.0)))) }).collect();

    let mut js = String::from("{\n");
    write!(js, "\"screen1\":\"{}\",\"screen2\":\"{}\",\"n_perm\":{nperm},\"n_boot\":{nboot},\"threads\":{},\"runs\":[\n", a[1], a[2], rayon::current_num_threads())?;
    let mut dens = String::from("umi_floor\tmatrix\tset\tpart\tx_lo\tx_hi\ty_lo\ty_hi\tcount\n");
    let mut first = true;
    for &floor in &[0.25f64, 0.5] {
        let (ix1, ix2) = shared_genes(&d1, &d2, floor);
        let n = ix1.len();
        let names: Vec<&str> = ix1.iter().map(|&i| d1.genes[i].as_str()).collect();
        let pos: HashMap<&str, usize> = names.iter().enumerate().map(|(i, g)| (*g, i)).collect();
        eprintln!("[{:.1}s] floor {floor}: N = {n}", t0.elapsed().as_secs_f64());
        let (f1, h0, h1, f2): (Mat, Mat, Mat, Mat) = (build(&d1, &d1.mean, &ix1, 1.0), build(&d1, &d1.h0, &ix1, 0.5), build(&d1, &d1.h1, &ix1, 0.5), build(&d2, &d2.mean, &ix2, 1.0));

        // pair classes
        let mut mask = vec![0u8; n * n];
        let mut tfpairs: Vec<(usize, usize)> = vec![]; // (TF, target), annotated in one direction only
        for i in 0..n {
            for j in 0..n {
                if i == j { continue; }
                if h0.t[i * n + j].abs().max(h0.t[j * n + i].abs()) >= 3.0 { mask[i * n + j] |= ACTIVE; }
                let fwd = tf.contains(&(names[i].to_string(), names[j].to_string()));
                let rev = tf.contains(&(names[j].to_string(), names[i].to_string()));
                if fwd || rev { mask[i * n + j] |= TF; }
                if fwd && !rev { tfpairs.push((i, j)); }
                if let (Some(x), Some(y)) = (tss.get(names[i]), tss.get(names[j])) {
                    if x.0 == y.0 && (x.1 - y.1).abs() < 1e6 { mask[i * n + j] |= CIS; }
                }
            }
        }
        for c in &complexes {
            let m: Vec<usize> = c.iter().filter_map(|g| pos.get(g.as_str()).copied()).collect();
            for &u in &m { for &v in &m { if u != v { mask[u * n + v] |= CPX; } } }
        }
        let id: Vec<usize> = (0..n).collect();
        for kind in ["scaled", "raw"] {
            let pick = |m: &Mat| if kind == "scaled" { m.sc.clone() } else { m.raw.clone() };
            let (mf1, mh1, mf2) = (pick(&f1), pick(&h1), pick(&f2));
            let ((sf1, af1), (sh1, ah1), (s2, a2)) = (split(&mf1, n), split(&mh1, n), split(&mf2, n));
            // label, held-out half of screen 1, required classes, excluded classes, permutation test
            let sets: [(&str, bool, u8, u8, bool); 9] = [
                ("all_pairs", false, 0, 0, true),
                ("active_pairs", true, ACTIVE, 0, true),
                ("all_pairs_without_neighbours", false, 0, CIS, false),
                ("active_pairs_without_neighbours", true, ACTIVE, CIS, false),
                ("active_pairs_outside_complexes", true, ACTIVE, CPX, false),
                ("tf_target_all", false, TF, 0, false),
                ("tf_target_active", true, TF | ACTIVE, 0, false),
                ("complex_all", false, CPX, 0, false),
                ("complex_active", true, CPX | ACTIVE, 0, false)];
            for (label, half, need, forbid, do_perm) in sets {
                let pr = Pairs { n, s1: if half { &sh1 } else { &sf1 }, a1: if half { &ah1 } else { &af1 }, s2: &s2, a2: &a2, mask: &mask, need, forbid };
                let obs = pr.stat(&id, &id);
                let (mut xs, mut ys, mut xa, mut ya) = (vec![], vec![], vec![], vec![]);
                for i in 0..n { for j in i + 1..n { if pr.keep(i, j) {
                    xs.push(pr.s1[i * n + j]); ys.push(pr.s2[i * n + j]); xa.push(pr.a1[i * n + j]); ya.push(pr.a2[i * n + j]); } } }
                let (rho_s, rho_a) = if xs.len() > 2 { (spearman(&xs, &ys), spearman(&xa, &ya)) } else { (f64::NAN, f64::NAN) };
                // binned scatter for the two headline sets
                if label == "all_pairs" || label == "active_pairs" {
                    for (part, x, y) in [("symmetric", &xs, &ys), ("skew", &xa, &ya)] {
                        let mut ab: Vec<f64> = x.iter().chain(y.iter()).map(|v| v.abs()).collect();
                        let lim = quant(&mut ab, 0.999).max(1e-9);
                        let nb = 60usize;
                        let mut h = vec![0u32; nb * nb];
                        let bin = |v: f64| (((v + lim) / (2.0 * lim) * nb as f64).floor().max(0.0) as usize).min(nb - 1);
                        for (u, v) in x.iter().zip(y.iter()) { h[bin(*u) * nb + bin(*v)] += 1; }
                        let w = 2.0 * lim / nb as f64;
                        for bx in 0..nb { for by in 0..nb { if h[bx * nb + by] > 0 {
                            writeln!(dens, "{floor}\t{kind}\t{label}\t{part}\t{:.5}\t{:.5}\t{:.5}\t{:.5}\t{}", -lim + bx as f64 * w, -lim + (bx + 1) as f64 * w, -lim + by as f64 * w, -lim + (by + 1) as f64 * w, h[bx * nb + by])?; } } }
                    }
                }
                let boots: Vec<Stat> = (0..nboot).into_par_iter().map(|b| {
                    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xB007 + b as u64);
                    let mut s: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
                    s.sort_unstable();
                    pr.stat(&s, &s)
                }).collect();
                let ci = |f: fn(&Stat) -> f64| { let mut v: Vec<f64> = boots.iter().map(f).collect(); let mut w = v.clone(); (quant(&mut v, 0.025), quant(&mut w, 0.975)) };
                let (ci_ra, ci_rs, ci_f) = (ci(|s| s.r_a), ci(|s| s.r_s), ci(|s| s.f_asym));
                let (mut p_ra, mut p_va, mut p_rs, mut null_ra_max) = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
                if do_perm && nperm > 0 {
                    let perms: Vec<Stat> = (0..nperm).into_par_iter().map(|b| {
                        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x9E37 + b as u64);
                        let mut p = id.clone();
                        p.shuffle(&mut rng);
                        pr.stat(&id, &p)
                    }).collect();
                    let pv = |f: fn(&Stat) -> f64| (1.0 + perms.iter().filter(|s| f(s) >= f(&obs)).count() as f64) / (1.0 + nperm as f64);
                    p_ra = pv(|s| s.r_a);
                    p_va = pv(|s| s.v_asym);
                    p_rs = pv(|s| s.r_s);
                    null_ra_max = perms.iter().map(|s| s.r_a).fold(f64::MIN, f64::max);
                }
                if !first { js.push_str(",\n"); }
                first = false;
                write!(js, "{{\"umi_floor\":{floor},\"matrix\":\"{kind}\",\"set\":\"{label}\",\"n_genes\":{n},\"n_pairs\":{},\"r_S\":{:.4},\"r_S_ci\":[{:.4},{:.4}],\"r_A\":{:.4},\"r_A_ci\":[{:.4},{:.4}],\"rho_S\":{:.4},\"rho_A\":{:.4},\"V_sym\":{:.4},\"V_asym\":{:.4},\"f_asym\":{:.4},\"f_asym_ci\":[{:.4},{:.4}],\"perm_p_r_A\":\"{:e}\",\"perm_p_V_asym\":\"{:e}\",\"perm_p_r_S\":\"{:e}\",\"perm_null_max_r_A\":\"{:.4}\"}}",
                    obs.npair, obs.r_s, ci_rs.0, ci_rs.1, obs.r_a, ci_ra.0, ci_ra.1, rho_s, rho_a, obs.v_sym, obs.v_asym, obs.f_asym, ci_f.0, ci_f.1, p_ra, p_va, p_rs, null_ra_max)?;
                eprintln!("[{:.1}s] {floor} {kind} {label}: pairs {} r_S {:.3} r_A {:.3} f_asym {:.3}", t0.elapsed().as_secs_f64(), obs.npair, obs.r_s, obs.r_a, obs.f_asym);
            }
            // direction of curated TF -> target pairs: select in one screen, read the direction in the other
            for (sel, mev, lab) in [(&f1, &mf2, "select_screen1_read_screen2"), (&f2, &mf1, "select_screen2_read_screen1")] {
                let (mut k_all, mut n_all, mut k_act, mut n_act) = (0usize, 0usize, 0usize, 0usize);
                for &(i, j) in &tfpairs {
                    let fwd = mev[i * n + j].abs() > mev[j * n + i].abs();
                    n_all += 1;
                    k_all += fwd as usize;
                    if sel.t[i * n + j].abs().max(sel.t[j * n + i].abs()) >= 3.0 { n_act += 1; k_act += fwd as usize; }
                }
                write!(js, ",\n{{\"umi_floor\":{floor},\"matrix\":\"{kind}\",\"set\":\"tf_direction_{lab}\",\"n_genes\":{n},\"n_tf_pairs\":{n_all},\"concordant_all\":{k_all},\"frac_all\":{:.4},\"p_all\":\"{:e}\",\"n_active\":{n_act},\"concordant_active\":{k_act},\"frac_active\":{:.4},\"p_active\":\"{:e}\"}}",
                    k_all as f64 / n_all.max(1) as f64, binom_p(k_all, n_all), k_act as f64 / n_act.max(1) as f64, binom_p(k_act, n_act))?;
            }
        }
    }
    write!(js, "\n],\"wall_seconds\":{:.1}\n}}\n", t0.elapsed().as_secs_f64())?;
    std::fs::write(format!("{}.json", a[6]), js)?;
    std::fs::write(format!("{}.density.tsv", a[6]), dens)?;
    eprintln!("done in {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}
