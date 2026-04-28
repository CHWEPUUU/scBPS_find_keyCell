#!/bin/bash
# set -euo pipefail

source /home/data/cq/.local/miniconda3/etc/profile.d/conda.sh 
conda activate scBPS 

gwas=${1:-Dutch207}
zscore_file=/home/data/share/module/scRNA/scBPS/GutMicrobiome_zscores/gene_zscore_${gwas}.txt

if [[ ! -f "${zscore_file}" ]]; then
    echo "[ERROR] zscore file not found: ${zscore_file}" >&2
    exit 1
fi

snakemake -s /home/data/share/module/scRNA/scBPS/code_module/scBPS.smk \
    --cores 32 \
    --config zscore_file="${zscore_file}" \
             scfile=./output/sc_adata.h5ad \
             outdir=./output/ \
             anno=./output/cell.annotation.txt