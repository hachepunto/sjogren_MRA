# Master Regulator Analysis of Primary Sjögren's Syndrome

Transcriptomic analysis of primary Sjögren's syndrome (pSS) in peripheral whole blood using Master Regulator Analysis (MRA) and differential expression across two independent GEO cohorts, with cross-dataset integration (Stouffer meta-analysis and Robust Rank Aggregation) and two supplementary validation cohorts.

---

## Cohorts

| | Discovery | Validation |
|---|---|---|
| **GEO accession** | [GSE51092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092) | [GSE84844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844) |
| **Platform** | Illumina HumanHT-12 v3.0 | Affymetrix HGU133 Plus 2.0 |
| **Design** | 32 controls + 190 pSS (222 total) | 30 controls + 30 pSS (60 total) |
| **Network method** | ARACNE-AP (server, 75 cores) | corto (Pearson, local) |

Supplementary validation: [GSE66795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66795) (29 controls + 131 pSS, Illumina HumanHT-12 v4.0; replication of DE genes) and [GSE140161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140161) (351 pSS, no controls; ssGSEA of the replicated up-regulated signature).

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
- `data/GSE66795_series_matrix.txt` (~87 MB)
- `data/GSE140161_series_matrix.txt` (~46 MB)

`ssgsea_GSE140161.R` additionally expects the processed gene-level matrix `data/GSE140161_Processed_data.txt.gz` from the GEO supplementary files.

### 2. Network inference (ARACNE-AP — server step)

The transcriptional network is inferred from pSS cases only. The full order:

1. Run **`mra_GSE84844.R`** to generate `results/networks/TFs_Lambert2018.txt`
2. Run **`mra_GSE51092.R`** up to step 3 — this exports `data/SJS_matrix.txt` (the cases-only expression matrix required by ARACNE-AP)
3. Set the `ARACNE` path in `run_aracne_ap.sh` to your [ARACNe-AP](https://github.com/califano-lab/ARACNe-AP) jar and run:
   ```bash
   bash run_aracne_ap.sh
   ```
   This requires ~300 GB RAM and 75 cores; adjust `N_CORES` and `MEM` as needed. Bootstraps go to `results/networks/aracne_ap_output/`; the consolidated network is copied to `results/networks/sjs_network.txt`.
4. Continue running **`mra_GSE51092.R`** from step 4 onward

### 3. R dependencies

```r
install.packages(c("BiocManager", "readr", "tibble", "purrr", "ggplot2", "patchwork",
                   "RobustRankAggreg", "openxlsx"))
BiocManager::install(c(
    "GEOquery", "matrixStats", "limma", "viper", "corto",
    "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "GSVA"
))
```

---

## Pipeline

Run scripts in order from the repository root (`Rscript script.R`, or source them from an R session whose working directory is the repo root). Every script creates the `results/` subfolders it needs, so no manual setup is required.

| Step | Script | Reads | Writes |
|------|--------|-------|--------|
| 0 | `get_data.sh` | — | `data/*_series_matrix.txt` |
| 1 | `mra_GSE84844.R` | `data/GSE84844_series_matrix.txt`, Lambert CSV | `results/networks/TFs_Lambert2018.txt`, `results/networks/regulon_GSE84844_corto.rda`, `results/mra/mra_GSE84844.rda`, `results/mra/top100_MRs_GSE84844.txt`, `results/de/DE_limma_GSE84844.txt`, `package_versions.tsv` |
| 2a | `mra_GSE51092.R` (steps 1–3) | `data/GSE51092_series_matrix.txt` | `data/SJS_matrix.txt` |
| 2b | `run_aracne_ap.sh` | `data/SJS_matrix.txt`, `results/networks/TFs_Lambert2018.txt` | `results/networks/sjs_network.txt` (server) |
| 2c | `mra_GSE51092.R` (steps 4–8) | `results/networks/sjs_network.txt` | `results/mra/mra_GSE51092.rda`, `results/mra/top100_MRs_GSE51092.txt`, `results/de/DE_limma_GSE51092.txt`, `results/figures/Top10_MRs_GSE51092.pdf` |
| 3 | `compare_MRs.R` | `results/mra/mra_*.rda` | `results/integration/MRs_comparacion.txt`, `results/integration/MRs_RRA_{activated,repressed}.txt`, `results/integration/comparison_results.rda`, `results/figures/MRs_scatter_NES.pdf` |
| 4 | `meta_MRs.R` | `results/integration/MRs_comparacion.txt`, `MRs_RRA_*.txt` | `results/integration/MRs_meta_shared.txt`, `results/integration/MRs_concordance_tests.txt`, `results/integration/meta_MRs_results.rda` |
| 5 | `figure_MRs.R` | `results/integration/comparison_results.rda`, `MRs_meta_shared.txt` | `results/figures/MRs_integrative_figure.pdf` |
| 6 | `ora_regulons.R` | `results/mra/mra_GSE51092.rda`, `results/integration/comparison_results.rda` | `results/enrichment/ORA_{GO_BP,Reactome}.txt`, `results/enrichment/ORA_results.rda`, `results/figures/ORA_combined.pdf` |
| 7 | `compare_ED.R` | `results/de/DE_limma_GSE{51092,84844}.txt` | `results/de/DE_comparacion.txt`, `results/de/DE_RRA_{up,down}regulated.txt`, `results/de/DE_replicated_{up,down}.txt`, `results/de/comparison_DE_results.rda`, `results/enrichment/DE_ORA_{GO_BP,Reactome}.txt`, `results/figures/DE_integrative_figure.pdf` |
| 8 | `meta_DE.R` | `results/de/DE_comparacion.txt`, `DE_RRA_*.txt` | `results/integration/DE_meta_shared.txt`, `results/integration/meta_DE_results.rda` |
| 9 | `convergence_MR_DE.R` | `results/integration/MRs_meta_shared.txt`, `DE_meta_shared.txt` | `results/integration/MR_DE_convergence.txt`, `results/integration/MR_DE_convergence_summary.txt` |
| 10 | `de_GSE66795.R` | `data/GSE66795_series_matrix.txt`, `results/de/comparison_DE_results.rda` | `results/validation/DE_limma_GSE66795.txt`, `results/validation/DE_GSE66795_replication.rda`, `results/validation/Table_S1_DE_GSE66795_replication.xlsx` |
| 11 | `ssgsea_GSE140161.R` | `data/GSE140161_Processed_data.txt.gz`, `results/de/DE_replicated_up.txt` | `results/validation/ssGSEA_GSE140161_scores_replicated.txt`, `results/validation/ssGSEA_GSE140161_replicated.rda`, `results/figures/ssGSEA_GSE140161_replicated.pdf` |
| 12 | `viper_cross_network.R` | `results/mra/mra_*.rda`, `results/networks/regulon_GSE84844_corto.rda`, `results/integration/MRs_meta_shared.txt` | `results/robustness/viper_cross_network{.txt,_concordance.txt,.rda}` |
| 13 | `viper_collectri.R` | `data/collectri_regulons_human.tsv`, `results/mra/mra_*.rda`, `MRs_meta_shared.txt` | `results/robustness/viper_collectri{.txt,_concordance.txt,.rda}`, `results/robustness/consensus_vs_collectri.txt` |
| 14 | `deconvolution_mra.R` | `results/mra/mra_*.rda`, corto regulon, `MRs_meta_shared.txt` | `results/robustness/xcell_*.txt`, `mcpcounter_*.txt`, `deconv_composition_tests.txt`, `mra_adjusted_composition.txt`, `mra_adjusted_concordance.txt`, `deconvolution_mra.rda` |
| 15 | `fgsea_regulons.R` | `results/mra/mra_*.rda`, corto regulon, `MRs_meta_shared.txt` | `results/enrichment/fgsea_regulons_{long,wide}.txt`, `fgsea_vs_viper_concordance.txt`, `fgsea_regulons.rda` |

Steps 3–5 form a chain: `compare_MRs.R` writes the shared-TF table that `meta_MRs.R` integrates, and `figure_MRs.R` draws the final figure from both outputs. Steps 7–8 are likewise sequential (`meta_DE.R` reads `DE_comparacion.txt` from `compare_ED.R`; `compare_ED.R` does not read anything from `meta_DE.R`). Step 9 needs both chains. Steps 6, 10 and 11 only need the outputs listed in their "Reads" column.

Note for OneDrive-synced checkouts: `pdf()` can occasionally fail at `dev.off()` with `write failed` while the sync client holds the file; rerun the script or write the figure to a local folder and copy it into `results/figures/`.

---

Steps 12–15 are robustness analyses added during peer review; they only need the `.rda` objects from steps 1–4. `data/collectri_regulons_human.tsv` is the CollecTRI TF–target table from OmniPath (`https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=1`, academic license; complexes removed, `mor = -1` for inhibition-only edges), and steps 13–14 additionally require `decoupleR`/`OmnipathR` (optional, only to re-download CollecTRI), `fgsea`, `xCell` (GitHub `dviraran/xCell`) and `MCPcounter` (GitHub `ebecht/MCPcounter`).

## Output layout (`results/`, git-ignored)

| Folder | Contents |
|--------|----------|
| `results/networks/` | TF list used by ARACNE-AP, consolidated ARACNE-AP network (`sjs_network.txt`), corto regulon (`regulon_GSE84844_corto.rda`) |
| `results/mra/` | `mra_GSE51092.rda`, `mra_GSE84844.rda`, `top100_MRs_*.txt` |
| `results/de/` | Per-cohort limma tables, shared-gene comparison, RRA lists, replicated gene lists, `comparison_DE_results.rda` |
| `results/integration/` | Cross-cohort MR integration (`MRs_meta_shared.txt`, `MRs_concordance_tests.txt`, `comparison_results.rda`, `meta_MRs_results.rda`), DE integration (`DE_meta_shared.txt`, `meta_DE_results.rda`), MR/DE convergence tables |
| `results/enrichment/` | ORA tables and objects for MR regulons (`ORA_*`) and for replicated DE genes (`DE_ORA_*`) |
| `results/validation/` | GSE66795 replication (limma table, `.rda`, Supplementary Table S1 `.xlsx`) and GSE140161 ssGSEA scores |
| `results/robustness/` | Cross-network VIPER, CollecTRI VIPER, composition-adjusted MRA (steps 12–14) |
| `results/figures/` | All PDF figures |

`package_versions.tsv` (written by the MRA scripts) stays in the repository root and is tracked for reproducibility.

### Key objects

| File | Contents |
|------|----------|
| `results/mra/mra_GSE51092.rda` | `mrs_gse51092`, `regulon_gse51092`, `eset`, `disease`, `de_gse51092`, `top_mrs` |
| `results/mra/mra_GSE84844.rda` | `mrs_gse84844`, `de_gse84844`, `top_mrs`, `eset`, `disease` (older files lack `eset`; `helpers_eset.R` rebuilds it on demand) |
| `results/integration/comparison_results.rda` | `shared` (TFs with a regulon in both networks), RRA results (`rra_act`, `rra_rep`) |
| `results/integration/meta_MRs_results.rda` | Stouffer + restricted RRA over shared TFs (`sh`, `concordance`) |
| `results/de/comparison_DE_results.rda` | `shared` genes with Stouffer z, replicated sets (`sig_up`, `sig_dn`), ORA objects |
| `results/integration/meta_DE_results.rda` | Stouffer + restricted RRA over shared genes (`d`) |
| `results/figures/MRs_integrative_figure.pdf` | Figure 2 — MR cross-dataset comparison |
| `results/figures/ORA_combined.pdf` | Figure 3 — ORA on MR target genes |
| `results/figures/DE_integrative_figure.pdf` | Figure 4 — DE cross-dataset + ORA |

---

## Citation

> [Manuscript in preparation]

---

## Data sources

- GSE51092: Becker et al. 2014 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092>
- GSE84844: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844>
- GSE66795: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66795>
- GSE140161: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140161>
- Lambert TF list: <https://humantfs.ccbr.utoronto.ca>
