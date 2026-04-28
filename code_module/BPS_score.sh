#!/bin/bash
# set -euo pipefail

gwas=${1:-Dutch207}
zscore_file=./GutMicrobiome_zscores/gene_zscore_${gwas}.txt

if [[ ! -f "${zscore_file}" ]]; then
    echo "[ERROR] zscore file not found: ${zscore_file}" >&2
    exit 1
fi

snakemake -s ./code_module/scBPS.smk \
    --cores 32 \
    --config zscore_file="${zscore_file}" \
             scfile=./output/sc_adata.h5ad \
             outdir=./output/ \
             anno=./output/cell.annotation.txt
