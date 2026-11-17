//! Read-only view of a sparse cell cache (genes.txt, groups.tsv, members.bin u32,
//! indptr.bin u64, lib.bin f32, idx.bin u16, cnt.bin u16, optional batch.bin u16).
use anyhow::{Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

pub struct Cache {
    pub genes: Vec<String>,
    pub groups: Vec<(String, usize, usize)>,
    members: Mmap,
    indptr: Mmap,
    idx: Mmap,
    cnt: Mmap,
    pub lib: Vec<f32>,
    pub batch: Vec<u16>,
}

fn map(p: &Path) -> Result<Mmap> {
    let f = File::open(p).with_context(|| format!("open {}", p.display()))?;
    Ok(unsafe { Mmap::map(&f)? })
}

impl Cache {
    pub fn open(dir: &str) -> Result<Self> {
        let d = Path::new(dir);
        let genes: Vec<String> = std::fs::read_to_string(d.join("genes.txt"))?.lines().map(String::from).collect();
        let groups = std::fs::read_to_string(d.join("groups.tsv"))?
            .lines()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())
            })
            .collect();
        let lib: Vec<f32> = bytemuck::cast_slice(&map(&d.join("lib.bin"))?).to_vec();
        let batch: Vec<u16> = match std::fs::read(d.join("batch.bin")) {
            Ok(b) => b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect(),
            Err(_) => vec![0; lib.len()],
        };
        Ok(Cache { genes, groups, members: map(&d.join("members.bin"))?, indptr: map(&d.join("indptr.bin"))?, idx: map(&d.join("idx.bin"))?, cnt: map(&d.join("cnt.bin"))?, lib, batch })
    }
    pub fn cells(&self, g: usize) -> &[u32] {
        let m: &[u32] = bytemuck::cast_slice(&self.members);
        let (_, o, n) = self.groups[g];
        &m[o..o + n]
    }
    #[inline]
    pub fn row(&self, c: usize) -> (&[u16], &[u16]) {
        let ip: &[u64] = bytemuck::cast_slice(&self.indptr);
        let (a, b) = (ip[c] as usize, ip[c + 1] as usize);
        let idx: &[u16] = bytemuck::cast_slice(&self.idx);
        let cnt: &[u16] = bytemuck::cast_slice(&self.cnt);
        (&idx[a..b], &cnt[a..b])
    }
}

pub fn read_f32(p: &str) -> Result<Vec<f32>> {
    let raw = std::fs::read(p).with_context(|| format!("read {p}"))?;
    Ok(raw.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
}
pub fn write_f32(p: &str, x: &[f32]) -> Result<()> {
    std::fs::write(p, bytemuck::cast_slice::<f32, u8>(x)).with_context(|| format!("write {p}"))
}

// ---------------------------------------------------------------------------------------
// Reduced screens: N x N mean and variance of standardised expression (row = knockdown).

use std::collections::HashMap;

pub struct Screen {
    pub genes: Vec<String>,
    pub ncell: Vec<f64>,
    pub umi: Vec<f64>,
    pub nctrl: f64,
    pub mean: Vec<f32>,
    pub var: Vec<f32>,
    /// two disjoint halves of the cells of each knockdown
    pub h0: Vec<f32>,
    pub h1: Vec<f32>,
}

fn read_sq(p: &str) -> Result<(Vec<String>, Vec<f64>, Vec<f64>, f64)> {
    let txt = std::fs::read_to_string(format!("{p}.sq.tsv")).with_context(|| format!("read {p}.sq.tsv"))?;
    let mut it = txt.lines();
    let nctrl: f64 = it.next().unwrap().split('\t').nth(1).unwrap().parse()?;
    let (mut genes, mut ncell, mut umi) = (vec![], vec![], vec![]);
    for l in it {
        let f: Vec<&str> = l.split('\t').collect();
        genes.push(f[0].to_string());
        ncell.push(f[1].parse()?);
        umi.push(f[2].parse()?);
    }
    Ok((genes, ncell, umi, nctrl))
}

impl Screen {
    pub fn load(p: &str) -> Result<Self> {
        let (genes, ncell, umi, nctrl) = read_sq(p)?;
        Ok(Screen { genes, ncell, umi, nctrl, mean: read_f32(&format!("{p}.mean.f32"))?, var: read_f32(&format!("{p}.var.f32"))?, h0: read_f32(&format!("{p}.h0.f32"))?, h1: read_f32(&format!("{p}.h1.f32"))? })
    }
    /// Two pseudo-screens from disjoint cell quarters of one screen: quarters 0 and 1, quarters 2 and 3.
    pub fn load_split(p: &str) -> Result<(Self, Self)> {
        let (genes, ncell, umi, nctrl) = read_sq(p)?;
        let q: Vec<Vec<f32>> = (0..4).map(|h| read_f32(&format!("{p}.q{h}.f32"))).collect::<Result<_>>()?;
        let avg = |x: &[f32], y: &[f32]| -> Vec<f32> { x.iter().zip(y).map(|(u, v)| 0.5 * (u + v)).collect() };
        let var = read_f32(&format!("{p}.var.f32"))?;
        let half: Vec<f64> = ncell.iter().map(|c| c / 2.0).collect();
        let a = Screen { genes: genes.clone(), ncell: half.clone(), umi: umi.clone(), nctrl, mean: avg(&q[0], &q[1]), var: var.clone(), h0: q[0].clone(), h1: q[1].clone() };
        let b = Screen { genes, ncell: half, umi, nctrl, mean: avg(&q[2], &q[3]), var, h0: q[2].clone(), h1: q[3].clone() };
        Ok((a, b))
    }
    /// The pair of screens named on the command line; "SPLIT" as second name means cell quarters of the first.
    pub fn pair(p1: &str, p2: &str) -> Result<(Self, Self)> {
        if p2 == "SPLIT" { Self::load_split(p1) } else { Ok((Self::load(p1)?, Self::load(p2)?)) }
    }
    pub fn n(&self) -> usize { self.genes.len() }
    /// t statistic of the shift of gene b under knockdown a in matrix m; frac = share of the cells behind m
    pub fn t(&self, m: &[f32], a: usize, b: usize, frac: f64) -> f64 {
        let n = self.n();
        let se = (self.var[a * n + b] as f64 / (self.ncell[a] * frac) + 1.0 / self.nctrl).sqrt();
        m[a * n + b] as f64 / se
    }
    pub fn se(&self, a: usize, b: usize) -> f64 {
        let n = self.n();
        (self.var[a * n + b] as f64 / self.ncell[a] + 1.0 / self.nctrl).sqrt()
    }
}

/// Genes knocked down with t <= -3 and expressed above the UMI floor in both screens.
pub fn shared_genes(d1: &Screen, d2: &Screen, floor: f64) -> (Vec<usize>, Vec<usize>) {
    let p2: HashMap<&str, usize> = d2.genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
    let (mut ix1, mut ix2) = (vec![], vec![]);
    for (i, g) in d1.genes.iter().enumerate() {
        let Some(&j) = p2.get(g.as_str()) else { continue };
        if d1.umi[i] < floor || d2.umi[j] < floor { continue; }
        if d1.t(&d1.mean, i, i, 1.0) > -3.0 || d2.t(&d2.mean, j, j, 1.0) > -3.0 { continue; }
        ix1.push(i);
        ix2.push(j);
    }
    (ix1, ix2)
}

/// N x N matrices for the kept genes: raw shift, knockdown-scaled response, t statistic.
pub struct Mat { pub raw: Vec<f64>, pub sc: Vec<f64>, pub t: Vec<f64> }
pub fn build(d: &Screen, m: &[f32], ix: &[usize], frac: f64) -> Mat {
    let (n, nn) = (ix.len(), d.n());
    let mut o = Mat { raw: vec![0.0; n * n], sc: vec![0.0; n * n], t: vec![0.0; n * n] };
    for (a, &ia) in ix.iter().enumerate() {
        let taa = d.t(m, ia, ia, frac);
        let w = taa * taa / (taa * taa + 1.0);
        let den = (m[ia * nn + ia] as f64).abs().max(1e-6);
        for (b, &ib) in ix.iter().enumerate() {
            let v = m[ia * nn + ib] as f64;
            o.raw[a * n + b] = v;
            o.sc[a * n + b] = v / den * w;
            o.t[a * n + b] = d.t(m, ia, ib, frac);
        }
    }
    o
}
/// Symmetric and skew-symmetric parts, diagonal set to zero.
pub fn split(m: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut s, mut a) = (vec![0f64; n * n], vec![0f64; n * n]);
    for i in 0..n { for j in 0..n { if i != j { s[i * n + j] = 0.5 * (m[i * n + j] + m[j * n + i]); a[i * n + j] = 0.5 * (m[i * n + j] - m[j * n + i]); } } }
    (s, a)
}
pub fn rank(x: &[f64]) -> Vec<f64> {
    let mut o: Vec<usize> = (0..x.len()).collect();
    o.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap());
    let mut r = vec![0f64; x.len()];
    let mut i = 0;
    while i < o.len() { let mut j = i; while j + 1 < o.len() && x[o[j + 1]] == x[o[i]] { j += 1; } let v = (i + j) as f64 / 2.0; for k in i..=j { r[o[k]] = v; } i = j + 1; }
    r
}
pub fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0f64, 0f64, 0f64);
    for (a, b) in x.iter().zip(y) { sxy += (a - mx) * (b - my); sxx += (a - mx).powi(2); syy += (b - my).powi(2); }
    sxy / (sxx * syy).sqrt()
}
pub fn spearman(x: &[f64], y: &[f64]) -> f64 { pearson(&rank(x), &rank(y)) }
/// exact two-sided binomial P against 0.5
pub fn binom_p(k: usize, n: usize) -> f64 {
    if n == 0 { return f64::NAN; }
    let mut lp = vec![0f64; n + 1];
    let mut c = 0f64;
    for i in 0..=n { if i > 0 { c += ((n - i + 1) as f64).ln() - (i as f64).ln(); } lp[i] = c - n as f64 * 2f64.ln(); }
    let obs = lp[k];
    lp.iter().filter(|&&l| l <= obs + 1e-9).map(|l| l.exp()).sum::<f64>().min(1.0)
}
