#!/usr/bin/env bash
# All analyses of the paper, from the reduced screens in reduced/ to the tables in results/.
# Run from the repository root after scripts/download_annotations.sh.
set -euo pipefail
cargo build --release
B=target/release
A=data/annotations
R=results
mkdir -p $R
bash scripts/complex_families.sh
LINES="k562ess rpe1 hepg2 jurkat"
rec() { $B/reciprocity "$1" "$2" $A/tf_dorothea_trrust.tsv $A/corum_complexes.tsv $A/tss.tsv "$3" 10000 2000 2> "$3.log"; }

# two independent K562 screens (genome-wide as discovery, essential-gene as replication)
rec reduced/gwps reduced/k562ess $R/k562_two_screens
$B/reciprocity reduced/gwps reduced/k562ess $A/tf_collectri.tsv $A/complex_families.tsv $A/tss.tsv $R/k562_two_screens_alt 0 2000 2> $R/k562_two_screens_alt.log
$B/calibration reduced/gwps reduced/k562ess 0.25 20 > $R/calibration.tsv
$B/gene_level reduced/gwps reduced/k562ess 0.25 $R/k562_two_screens.gene_level
$B/directed_pairs reduced/gwps reduced/k562ess 0.25 $R/k562_two_screens.directed
$B/linear_response reduced/k562ess reduced/gwps 0.25 > $R/linear_response.tsv

# within each cell line: disjoint cell quarters of one screen
for l in $LINES; do
  rec reduced/$l SPLIT $R/within_$l
  $B/reciprocity reduced/$l SPLIT $A/tf_collectri.tsv $A/complex_families.tsv $A/tss.tsv $R/within_${l}_alt 0 2000 2> $R/within_${l}_alt.log
  $B/gene_level reduced/$l SPLIT 0.25 $R/within_$l.gene_level
  $B/directed_pairs reduced/$l SPLIT 0.25 $R/within_$l.directed
  $B/linear_response reduced/$l SPLIT 0.25 | tail -n +2 >> $R/linear_response.tsv
done

# between cell lines
set -- $LINES
while [ $# -gt 1 ]; do
  a=$1; shift
  for b in "$@"; do rec reduced/$a reduced/$b $R/between_${a}_$b; done
done
echo "all analyses finished"
