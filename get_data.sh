#!/usr/bin/env bash
# Downloads and decompresses GEO expression matrices required by this pipeline.
# Run once from the repo root before executing any R scripts.
set -euo pipefail

DATA_DIR="$(cd "$(dirname "$0")/data" && pwd)"

GSE51092_GZ="$DATA_DIR/GSE51092_series_matrix.txt.gz"
GSE84844_GZ="$DATA_DIR/GSE84844_series_matrix.txt.gz"

if [[ ! -f "$DATA_DIR/GSE51092_series_matrix.txt" ]]; then
    curl -L --progress-bar \
        -o "$GSE51092_GZ" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51092/matrix/GSE51092_series_matrix.txt.gz"
    gunzip "$GSE51092_GZ"
fi

if [[ ! -f "$DATA_DIR/GSE84844_series_matrix.txt" ]]; then
    curl -L --progress-bar \
        -o "$GSE84844_GZ" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84844/matrix/GSE84844_series_matrix.txt.gz"
    gunzip "$GSE84844_GZ"
fi
