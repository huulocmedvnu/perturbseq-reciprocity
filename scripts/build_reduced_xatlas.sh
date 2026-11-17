#!/usr/bin/env bash
# X-Atlas/Orion genome-wide screens (HCT116, HEK293T) from the public parquet shards to reduced/.
# Run from the repository root. Needs about 130 GB of downloads (deleted shard by shard), 60 GB of
# disk for the caches, 64 GB of memory, Rust, curl and Python 3.
#   Data: Huang et al. 2025, https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion (CC BY-NC-SA 4.0)
# One shard is one GEM batch. Cells that pass the guide filter are kept. Genes are kept when they are
# knockdown targets with a mean of at least 0.25 UMI per control cell, the lowest floor of the analyses.
#   JOBS            shards downloaded and converted at the same time (default 8)
#   REDUCE_THREADS  threads of the reduce step, which holds one gene-by-gene matrix per thread (default 8)
set -euo pipefail
H=https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/main
API="https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Orion/tree/main/data?limit=1000"
W=raw/xatlas
mkdir -p $W cache reduced
cargo build --release
(cd tools/xatlas_cache && cargo build --release)
export X=$PWD/tools/xatlas_cache/target/release/xatlas_cache H W
[ -s $W/gene_metadata.parquet ] || curl -sfL --retry 20 -o $W/gene_metadata.parquet $H/metadata/gene_metadata.parquet
[ -s $W/shards.txt ] || curl -sf --retry 20 "$API" | python3 -I -c 'import json, sys
for e in json.load(sys.stdin): print(e["path"].split("/")[-1])' > $W/shards.txt
one() {
  local s=$1 p=$W/${1%.parquet}
  [ -f $p.done ] && return 0
  curl -sfL -C - --retry 20 --retry-delay 5 --speed-limit 100000 --speed-time 60 -o $p.parquet $H/data/$s && $X shard $p.parquet $p && rm -f $p.parquet
}
export -f one
for l in HCT116 HEK293T; do
  o=$(echo $l | tr A-Z a-z)
  [ -f reduced/$o.ctrlcorr.f32 ] && continue
  grep "^${l}_" $W/shards.txt > $W/$l.txt
  xargs -P "${JOBS:-8}" -I{} bash -c 'one {}' < $W/$l.txt
  $X merge $W $W/$l.txt $W/gene_metadata.parquet 0.25 cache/$o
  RAYON_NUM_THREADS="${REDUCE_THREADS:-8}" target/release/reduce cache/$o reduced/$o
  rm -f $W/${l}_Batch*
done
ls -la reduced
