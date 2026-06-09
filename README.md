# Master Regulator Analysis of Primary Sjögren's Syndrome

Transcriptomic analysis of primary Sjögren's syndrome (pSS) in peripheral whole blood
using Master Regulator Analysis (MRA) and differential expression across two independent
GEO cohorts, with cross-dataset validation by Robust Rank Aggregation.

---

## Cohorts

| | Discovery | Validation |
|---|---|---|
| **GEO accession** | [GSE51092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092) | [GSE84844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844) |
| **Platform** | Illumina HumanHT-12 v3.0 | Affymetrix HGU133 Plus 2.0 |
| **Design** | 32 controls + 190 pSS (222 total) | 30 controls + 30 pSS (60 total) |
| **Network method** | ARACNE-AP (server, 75 cores) | corto (Pearson, local) |

Transcription factor list: Lambert et al. 2018 (Cell, PMID: 29425488) — 1,639 confirmed human TFs.

---

## Setup

### 1. Download expression matrices

```bash
bash get_data.sh
```

Downloads and decompresses the GEO series matrix files into `data/`:
- `data/GSE51092_series_matrix.txt` (~39 MB)
- `data/GSE84844_series_matrix.txt` (~38 MB)

### 2. Network inference (ARACNE-AP — server step)

The transcriptional network is inferred from pSS cases only. The full order:

1. Run **`mra_GSE84844.R`** to generate `data/TFs_Lambert2018.txt`
2. Run **`mra_GSE51092.R`** up to step 3 — this exports `data/SJS_matrix.txt`
   (the cases-only expression matrix required by ARACNE-AP)
3. Set the `ARACNE` path in `run_aracne_ap.sh` to your
   [ARACNe-AP](https://github.com/califano-lab/ARACNe-AP) jar and run:
   ```bash
   bash run_aracne_ap.sh
   ```
   This requires ~300 GB RAM and 75 cores; adjust `N_CORES` and `MEM` as needed.
   Output: `data/sjs_network.txt`
4. Continue running **`mra_GSE51092.R`** from step 4 onward

### 3. R dependencies

```r
install.packages(c("BiocManager", "readr", "tibble", "purrr", "ggplot2", "patchwork"))
BiocManager::install(c(
    "GEOquery", "matrixStats", "limma", "viper", "corto",
    "clusterProfiler", "ReactomePA", "org.Hs.eg.db"
))
install.packages("RobustRankAggreg")
```

---

## Pipeline

Run scripts in order from the repository root:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `mra_GSE84844.R` | Downloads Lambert TF list → `data/`; corto network + MRA, validation cohort |
| 2a | `mra_GSE51092.R` (steps 1–3) | Exports `data/SJS_matrix.txt` (cases-only matrix for ARACNE-AP) |
| 2b | `run_aracne_ap.sh` | ARACNE-AP network inference → `data/sjs_network.txt` (server) |
| 2c | `mra_GSE51092.R` (steps 4–9) | Loads network, runs MRA + DE, discovery cohort |
| 3 | `compare_MRs.R` | Cross-dataset MR comparison + RRA → integrative figure |
| 4 | `ora_regulons.R` | Over-representation analysis on MR target genes |
| 5 | `compare_ED.R` | Differential expression + RRA + ORA → integrative figure |

Each script is self-contained and writes its outputs (`.rda`, `.txt`, `.pdf`) to the
working directory. Set your working directory to the repository root before running,
or open `sjogren_MRA.Rproj` if using RStudio.

---

## Key outputs

| File | Contents |
|------|----------|
| `mra_GSE51092.rda` | `mrs_gse51092`, `regulon_gse51092`, `eset`, `disease`, `de_gse51092`, `top_mrs` |
| `mra_GSE84844.rda` | `mrs_gse84844`, `de_gse84844`, `top_mrs` |
| `comparison_results.rda` | RRA results for MRs (`rra_act`, `rra_rep`) |
| `comparison_DE_results.rda` | RRA + ORA results for DE genes |
| `MRs_integrative_figure.pdf` | Figure 2 — MR cross-dataset comparison |
| `ORA_combined.pdf` | Figure 3 — ORA on MR target genes |
| `DE_integrative_figure.pdf` | Figure 4 — DE cross-dataset + ORA |

---

## Citation

> [Manuscript in preparation]

---

## Data sources

- GSE51092: Becker et al. 2014 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092>
- GSE84844: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844>
- Lambert TF list: <https://humantfs.ccbr.utoronto.ca>
