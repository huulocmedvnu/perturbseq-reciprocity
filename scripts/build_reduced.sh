#!/usr/bin/env bash
# From the public single-cell files to the reduced screens in reduced/, and to the shifts of unnormalised
# counts in reduced_raw/ for the four essential-gene screens. Run from the repository root.
# Needs about 40 GB of disk, the HDF5 library, Rust, Python 3 with h5py.
#   K562 genome-wide, K562 essential, RPE1 : Replogle et al. 2022, as distributed by scPerturb (Zenodo 13350497)
#   HepG2, Jurkat                          : Nadig et al. 2025 (GEO GSE264667)
# The sparse cell cache is written by the cache step of dosecurve (https://github.com/huulocmedvnu/dosecurve, v1.0.0).
set -euo pipefail
mkdir -p raw cache reduced reduced_raw builder
[ -d builder/dosecurve ] || git clone --depth 1 --branch v1.0.0 https://github.com/huulocmedvnu/dosecurve builder/dosecurve
(cd builder/dosecurve/tools/dosecurve && cargo build --release)
DC=builder/dosecurve/tools/dosecurve/target/release/dosecurve
T=builder/dosecurve/tools
cargo build --release
Z=https://zenodo.org/records/13350497/files
G=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl
for pair in "gwps ReplogleWeissman2022_K562_gwps" "k562ess ReplogleWeissman2022_K562_essential" "rpe1 ReplogleWeissman2022_rpe1"; do
  set -- $pair
  [ -f reduced/$1.ctrlcorr.f32 ] && continue
  curl -L -C - --retry 20 -o raw/$1.h5ad "$Z/$2.h5ad?download=1"
  $DC cache raw/$1.h5ad cache/$1
  target/release/reduce cache/$1 reduced/$1
  if [ $1 != gwps ]; then target/release/reduce cache/$1 reduced_raw/$1 raw; fi
  rm -f raw/$1.h5ad
done
for l in hepg2 jurkat; do
  [ -f reduced/$l.ctrlcorr.f32 ] && continue
  curl -L -C - --retry 20 -o raw/$l.h5ad "$G/GSE264667_${l}_raw_singlecell_01.h5ad"
  python3 -I $T/h5ad_addcol.py raw/$l.h5ad
  $DC cache raw/$l.h5ad cache/$l
  python3 -I $T/genes_from_var.py raw/$l.h5ad cache/$l
  target/release/reduce cache/$l reduced/$l
  target/release/reduce cache/$l reduced_raw/$l raw
  rm -f raw/$l.h5ad
done
ls -la reduced
