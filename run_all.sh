#!/usr/bin/env bash
# All analyses of the paper, from the reduced screens in reduced/ to the tables in results/.
# Run from the repository root after scripts/download_annotations.sh.
set -euo pipefail
cargo build --release
B=target/release
A=data/annotations
R=results
mkdir -p $R
rm -f $R/transfer_prediction.tsv
bash scripts/complex_families.sh
LINES="k562ess rpe1 hepg2 jurkat hct116 hek293t"   # screens with cell quarters
ESS="k562ess rpe1 hepg2 jurkat"                    # essential-gene screens, used for the prediction analyses
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

# prediction of held-out knockdowns: standardised response with the other cell lines as reference,
# then the shift of unnormalised counts (reduced_raw/, written by scripts/build_reduced.sh)
P="0.25 $A/tf_dorothea_trrust.tsv $A/corum_complexes.tsv 5"
others() { local o=""; for r in $ESS; do [ $r != $1 ] && o="$o,reduced/$r"; done; echo "refs=${o#,}"; }
$B/heldout_prediction reduced/k562ess reduced/gwps $P $(others k562ess) > $R/heldout_prediction.tsv
for l in $ESS; do $B/heldout_prediction reduced/$l SPLIT $P $(others $l) | tail -n +2 >> $R/heldout_prediction.tsv; done
for l in $ESS; do $B/heldout_prediction reduced/$l SPLIT $P reduced_raw/$l | tail -n +2 >> $R/heldout_prediction.tsv; done
# response of one cell line from the screens of the other three
for t in $ESS; do
  refs=""; for r in $ESS; do [ $r != $t ] && refs="$refs reduced/$r"; done
  $B/transfer_prediction 0.25 reduced/$t $refs | { if [ $t = k562ess ]; then cat; else tail -n +2; fi; } >> $R/transfer_prediction.tsv
done

# between cell lines
set -- $LINES
while [ $# -gt 1 ]; do
  a=$1; shift
  for b in "$@"; do rec reduced/$a reduced/$b $R/between_${a}_$b; done
done
# genome-wide K562 screen against the genome-wide screens of two other cell lines
for b in hct116 hek293t; do rec reduced/gwps reduced/$b $R/between_gwps_$b; done
echo "all analyses finished"
