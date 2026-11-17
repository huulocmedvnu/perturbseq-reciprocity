#!/usr/bin/env bash
# Curated annotations used to label gene pairs. Run from the repository root.
#   transcription factor -> target : DoRothEA (confidence A and B) and TRRUST v2, CollecTRI
#   protein complexes               : CORUM, through OmniPath
#   transcription start sites       : RefSeq Select / MANE transcripts (UCSC hg38 ncbiRefSeqSelect)
set -euo pipefail
A=data/annotations
mkdir -p $A
curl -sL "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv" -o $A/trrust.raw.tsv
curl -sL "https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B&genesymbols=1&fields=dorothea_level&license=academic" -o $A/dorothea.raw.tsv
curl -sL "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&license=academic" -o $A/collectri.raw.tsv
curl -sL "https://omnipathdb.org/complexes?resources=CORUM&license=academic" -o $A/corum.raw.tsv
{ awk -F'\t' '{print $1"\t"$2}' $A/trrust.raw.tsv; awk -F'\t' 'NR>1{print $3"\t"$4}' $A/dorothea.raw.tsv; } | sort -u > $A/tf_dorothea_trrust.tsv
awk -F'\t' 'NR>1{print $3"\t"$4}' $A/collectri.raw.tsv | sort -u > $A/tf_collectri.tsv
awk -F'\t' 'NR>1 && $3!=""{print $1"\t"$3}' $A/corum.raw.tsv > $A/corum_complexes.tsv
curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz" | gunzip -c \
  | awk -F'\t' '{tss = ($4=="+") ? $5+1 : $6; print $13"\t"$3"\t"tss}' | sort -u -k1,1 > $A/tss.tsv
rm -f $A/*.raw.tsv
wc -l $A/*.tsv
