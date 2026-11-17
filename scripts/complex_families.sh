#!/usr/bin/env bash
# Six multiprotein complexes defined by gene-name family, from the genes of the reduced screens.
# Used as a second definition of complex subunits next to CORUM.
set -euo pipefail
cat reduced/*.sq.tsv | cut -f1 | grep -v '^#' | sort -u | awk '
  /^RP[LS][0-9]/ || /^RPLP[0-2]$/ || /^RPSA$/ {a["cytoplasmic ribosome"]=a["cytoplasmic ribosome"]"_"$1; next}
  /^MRP[LS][0-9]/ {a["mitochondrial ribosome"]=a["mitochondrial ribosome"]"_"$1; next}
  /^PSM[ABCD][0-9]/ {a["proteasome"]=a["proteasome"]"_"$1; next}
  /^POLR2/ {a["RNA polymerase II"]=a["RNA polymerase II"]"_"$1; next}
  /^MED[0-9]/ {a["Mediator"]=a["Mediator"]"_"$1; next}
  /^SF3[AB][0-9]/ || /^SNRP/ || /^LSM[0-9]/ || /^PRPF[0-9]/ {a["spliceosome"]=a["spliceosome"]"_"$1; next}
  END {for (k in a) print k"\t"substr(a[k],2)}' > data/annotations/complex_families.tsv
