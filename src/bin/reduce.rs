//! Reduce a sparse cell cache to the square response summary.
//! Expression x = ln(1 + CP10k). Every cell is standardised against the control cells of
//! its own batch: z = (x - mu_batch) / sd_batch. For each knocked-down gene A that is also
//! measured, the mean and variance of z of every square gene B across the cells of A are
//! written, for all cells, for two interleaved halves (h0, h1) and for four interleaved
//! quarters (q0..q3; h0 = q0 + q2, h1 = q1 + q3). The correlation matrix of z among the
//! square genes in control cells is written as well.
//! Usage: reduce <cache_dir> <out_prefix> [raw]
//! With "raw" the unnormalised counts are used without batches and without scaling.
//! Output: <out>.sq.tsv (gene, n_cells, mean control UMI per cell),
//!         <out>.mean.f32 <out>.var.f32 <out>.h0.f32 <out>.h1.f32 <out>.q0..q3.f32
//!         (N x N, row = knockdown), <out>.ctrlcorr.f32 (N x N)
use anyhow::{bail, Result};
use perturbseq_reciprocity::{write_f32, Cache};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

const SD_FLOOR: f64 = 0.05;
const NONE: u32 = u32::MAX;

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 3 && a.len() != 4 { bail!("usage: reduce <cache_dir> <out_prefix> [raw]"); }
    // raw: unnormalised counts, no batches, no scaling. The files then hold shifts of mean counts,
    // and <out>.ctrlsd.f32 holds the standard deviation of counts in control cells.
    let raw = a.get(3).map(|s| s == "raw").unwrap_or(false);
    let tx = move |v: f64, sc: f64| -> f64 { if raw { v } else { (1.0 + v * sc).ln() } };
    let t = Instant::now();
    let c = Cache::open(&a[1])?;
    let col: HashMap<&str, usize> = c.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
    let ctrl = c.groups.iter().position(|g| g.0 == "control").expect("control group");
    // square set: knocked-down genes that are also measured
    let sq: Vec<(usize, usize)> = c.groups.iter().enumerate()
        .filter_map(|(g, x)| col.get(x.0.as_str()).map(|&j| (g, j))).collect();
    let n = sq.len();
    let mut cmap = vec![NONE; c.genes.len()];
    for (k, &(_, j)) in sq.iter().enumerate() { cmap[j] = k as u32; }
    let nb = if raw { 1 } else { *c.batch.iter().max().unwrap() as usize + 1 };
    let bat = |cell: usize| -> usize { if raw { 0 } else { c.batch[cell] as usize } };
    eprintln!("{} cells, {} genes, {} groups, square N = {}, batches = {}", c.lib.len(), c.genes.len(), c.groups.len(), n, nb);

    // control statistics per batch
    let mut cnt_b = vec![0f64; nb];
    let mut s1 = vec![0f64; nb * n];
    let mut s2 = vec![0f64; nb * n];
    let mut umi = vec![0f64; n];
    for &cell in c.cells(ctrl) {
        let cell = cell as usize;
        let b = bat(cell);
        let sc = 1e4 / c.lib[cell] as f64;
        cnt_b[b] += 1.0;
        let (ix, ct) = c.row(cell);
        for (&j, &v) in ix.iter().zip(ct) {
            let k = cmap[j as usize];
            if k == NONE { continue; }
            let x = tx(v as f64, sc);
            s1[b * n + k as usize] += x;
            s2[b * n + k as usize] += x * x;
            umi[k as usize] += v as f64;
        }
    }
    let n_ctrl = c.cells(ctrl).len() as f64;
    // mu, 1/sd and z of a zero count, per batch and gene
    let mut mu = vec![0f32; nb * n];
    let mut isd = vec![0f32; nb * n];
    let mut z0 = vec![0f32; nb * n];
    let mut csd = vec![0f32; n];
    for b in 0..nb {
        if cnt_b[b] < 20.0 { continue; } // batch without usable controls: its cells get z = 0
        for k in 0..n {
            let m = s1[b * n + k] / cnt_b[b];
            let v = (s2[b * n + k] / cnt_b[b] - m * m).max(0.0) * cnt_b[b] / (cnt_b[b] - 1.0);
            if v <= 0.0 { continue; }
            let sd = if raw { csd[k] = v.sqrt() as f32; 1.0 } else { v.sqrt().max(SD_FLOOR) };
            mu[b * n + k] = m as f32;
            isd[b * n + k] = (1.0 / sd) as f32;
            z0[b * n + k] = (-m / sd) as f32;
        }
    }
    eprintln!("controls done in {:.1}s", t.elapsed().as_secs_f64());

    let mut mean = vec![0f32; n * n];
    let mut var = vec![0f32; n * n];
    let mut h0 = vec![0f32; n * n];
    let mut h1 = vec![0f32; n * n];
    let mut qq = vec![0f32; 4 * n * n]; // row-wise: for knockdown k, the four quarters one after the other
    (mean.par_chunks_mut(n), var.par_chunks_mut(n), h0.par_chunks_mut(n), h1.par_chunks_mut(n), qq.par_chunks_mut(4 * n), sq.par_iter())
        .into_par_iter()
        .for_each(|(rm, rv, r0, r1, rq, &(g, _))| {
            let mut a1 = [vec![0f64; n], vec![0f64; n], vec![0f64; n], vec![0f64; n]];
            let mut a2 = vec![0f64; n];
            let mut nbh = [vec![0f64; nb], vec![0f64; nb], vec![0f64; nb], vec![0f64; nb]];
            for (pos, &cell) in c.cells(g).iter().enumerate() {
                let cell = cell as usize;
                let b = bat(cell);
                let h = pos & 3;
                nbh[h][b] += 1.0;
                let sc = 1e4 / c.lib[cell] as f64;
                let (ix, ct) = c.row(cell);
                let (mub, isdb, z0b) = (&mu[b * n..(b + 1) * n], &isd[b * n..(b + 1) * n], &z0[b * n..(b + 1) * n]);
                for (&j, &v) in ix.iter().zip(ct) {
                    let k = cmap[j as usize];
                    if k == NONE { continue; }
                    let k = k as usize;
                    let z = (tx(v as f64, sc) - mub[k] as f64) * isdb[k] as f64;
                    let zz = z0b[k] as f64;
                    a1[h][k] += z - zz;
                    a2[k] += z * z - zz * zz;
                }
            }
            // add the contribution of zero counts: every cell of batch b carries z0_b
            let mut nq = [0f64; 4];
            for b in 0..nb {
                let c: [f64; 4] = [nbh[0][b], nbh[1][b], nbh[2][b], nbh[3][b]];
                let ct: f64 = c.iter().sum();
                if ct == 0.0 { continue; }
                let z0b = &z0[b * n..(b + 1) * n];
                for h in 0..4 { nq[h] += c[h]; }
                for k in 0..n {
                    let zz = z0b[k] as f64;
                    for h in 0..4 { a1[h][k] += c[h] * zz; }
                    a2[k] += ct * zz * zz;
                }
            }
            let nt: f64 = nq.iter().sum();
            let (n0, n1) = (nq[0] + nq[2], nq[1] + nq[3]);
            for k in 0..n {
                let m = (a1[0][k] + a1[1][k] + a1[2][k] + a1[3][k]) / nt;
                rm[k] = m as f32;
                rv[k] = if nt > 1.0 { ((a2[k] - nt * m * m) / (nt - 1.0)).max(0.0) as f32 } else { 0.0 };
                r0[k] = if n0 > 0.0 { ((a1[0][k] + a1[2][k]) / n0) as f32 } else { 0.0 };
                r1[k] = if n1 > 0.0 { ((a1[1][k] + a1[3][k]) / n1) as f32 } else { 0.0 };
                for h in 0..4 { rq[h * n + k] = if nq[h] > 0.0 { (a1[h][k] / nq[h]) as f32 } else { 0.0 }; }
            }
        });
    let o = &a[2];
    if raw { write_f32(&format!("{o}.ctrlsd.f32"), &csd)?; }
    write_f32(&format!("{o}.mean.f32"), &mean)?;
    write_f32(&format!("{o}.var.f32"), &var)?;
    write_f32(&format!("{o}.h0.f32"), &h0)?;
    write_f32(&format!("{o}.h1.f32"), &h1)?;
    for h in 0..4 {
        let mut m = vec![0f32; n * n];
        for k in 0..n { m[k * n..(k + 1) * n].copy_from_slice(&qq[k * 4 * n + h * n..k * 4 * n + (h + 1) * n]); }
        write_f32(&format!("{o}.q{h}.f32"), &m)?;
    }
    drop(qq);
    eprintln!("knockdowns done in {:.1}s", t.elapsed().as_secs_f64());

    // correlation of z among the square genes in control cells.
    // z = z0_b + d, with d non-zero only where the count is non-zero, so
    // sum z z' = sum_b [ n_b z0 z0' + z0 D_b' + D_b z0' ] + sum_cells d d'.
    let ctrl_cells = c.cells(ctrl);
    let (dd, dsum) = ctrl_cells.par_chunks(256).fold(
        || (vec![0f64; n * n], vec![0f64; nb * n]),
        |(mut acc, mut ds), chunk| {
            let mut kk: Vec<usize> = Vec::new();
            let mut dv: Vec<f64> = Vec::new();
            for &cell in chunk {
                let cell = cell as usize;
                let b = bat(cell);
                if cnt_b[b] < 20.0 { continue; }
                let sc = 1e4 / c.lib[cell] as f64;
                let (ix, ct) = c.row(cell);
                kk.clear(); dv.clear();
                for (&j, &v) in ix.iter().zip(ct) {
                    let k = cmap[j as usize];
                    if k == NONE { continue; }
                    let k = k as usize;
                    let d = (tx(v as f64, sc) - mu[b * n + k] as f64) * isd[b * n + k] as f64 - z0[b * n + k] as f64;
                    kk.push(k); dv.push(d);
                    ds[b * n + k] += d;
                }
                for (x, &kx) in kk.iter().enumerate() {
                    let row = &mut acc[kx * n..(kx + 1) * n];
                    let dx = dv[x];
                    for (y, &ky) in kk.iter().enumerate() { row[ky] += dx * dv[y]; }
                }
            }
            (acc, ds)
        }).reduce(
        || (vec![0f64; n * n], vec![0f64; nb * n]),
        |(mut a, mut b), (c2, d2)| { for (x, y) in a.iter_mut().zip(&c2) { *x += y; } for (x, y) in b.iter_mut().zip(&d2) { *x += y; } (a, b) });
    let mut cc = dd;
    cc.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        for b in 0..nb {
            if cnt_b[b] < 20.0 { continue; }
            let (zi, di) = (z0[b * n + i] as f64, dsum[b * n + i]);
            let (zb, db) = (&z0[b * n..(b + 1) * n], &dsum[b * n..(b + 1) * n]);
            let w = cnt_b[b] * zi;
            for j in 0..n { row[j] += w * zb[j] as f64 + zi * db[j] + di * zb[j] as f64; }
        }
    });
    let diag: Vec<f64> = (0..n).map(|i| cc[i * n + i]).collect();
    let corr: Vec<f32> = (0..n * n).map(|x| { let (i, j) = (x / n, x % n); let d = (diag[i] * diag[j]).sqrt(); if d > 0.0 { (cc[x] / d) as f32 } else { 0.0 } }).collect();
    write_f32(&format!("{o}.ctrlcorr.f32"), &corr)?;
    let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{o}.sq.tsv"))?);
    writeln!(f, "#n_control\t{}", n_ctrl)?;
    for (k, &(g, j)) in sq.iter().enumerate() {
        writeln!(f, "{}\t{}\t{:.5}", c.genes[j], c.groups[g].2, umi[k] / n_ctrl)?;
    }
    eprintln!("done in {:.1}s", t.elapsed().as_secs_f64());
    Ok(())
}
