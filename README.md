# Master Regulator Analysis of Primary Sjögren's Syndrome

Transcriptomic analysis of primary Sjögren's syndrome (pSS) in peripheral whole blood using Master Regulator Analysis (MRA) and differential expression across two independent GEO cohorts, with cross-dataset integration (Stouffer meta-analysis of regulator activity and replication of differential expression), robustness analyses (exchanged networks, CollecTRI, cell-composition adjustment, fgsea) and two supplementary validation cohorts. The whole pipeline is a Snakemake workflow.

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

### R environment

R 4.6 with Bioconductor 3.23. `envs/install_r_packages.R` installs every package the pipeline needs (CRAN, Bioconductor, corto 1.2.4 and the GitHub packages `xCell` and `MCPcounter`); `package_versions.tsv` lists the exact versions used for the manuscript. With conda, for example:

```bash
conda create -n sjogren -c conda-forge r-base=4.6.1 c-compiler cxx-compiler fortran-compiler make cmake pkg-config \
    libcurl libxml2 openssl zlib icu glpk gmp freetype fontconfig harfbuzz fribidi libpng libjpeg-turbo libtiff cairo imagemagick
conda activate sjogren      # packages must compile with the environment's own toolchain
Rscript envs/install_r_packages.R
```

The workflow itself needs [Snakemake](https://snakemake.readthedocs.io) (≥ 8) and `curl`; Java and the [ARACNe-AP](https://github.com/califano-lab/ARACNe-AP) jar only if the network is re-inferred.

### Configuration (`config.yaml`)

| Key | Default | Meaning |
|-----|---------|---------|
| `rscript` | `Rscript` | Path to the `Rscript` of the R 4.6 installation |
| `corto_threads` | `4` | Threads for corto (GSE84844 network) |
| `run_aracne` | `false` | `false`: use the frozen network in `resources/`; `true`: re-infer it with ARACNe-AP |
| `aracne_jar` | — | Path to `aracne.jar` (only with `run_aracne: true`) |
| `aracne_bootstraps`, `aracne_pvalue`, `aracne_mem` | `100`, `1E-8`, `4G` | ARACNe-AP parameters; each bootstrap is one job using `aracne_mem` of RAM |

---

## Running

```bash
snakemake --cores 8                         # full pipeline with the frozen ARACNe-AP network
snakemake --cores 75 --config run_aracne=true aracne_jar=/path/to/aracne.jar   # re-infer the network
snakemake -n                                # dry run: list pending jobs
```

Snakemake downloads the GEO series matrices, the GSE140161 processed matrix and the Lambert et al. TF table into `data/`, runs every step in dependency order and writes one log per step to `logs/`. Every script can also be run on its own from the repository root (`Rscript scripts/<name>.R`) once its inputs exist.

### Frozen resources (`resources/`)

| File | Contents |
|------|----------|
| `sjs_network_aracne.txt.gz` | Consolidated ARACNe-AP network of GSE51092 (190 pSS cases, Lambert TFs, p < 1e-8, 100 bootstraps with seeds 1–100) |
| `probe_map_GSE51092_2019.tsv` | Gene → Illumina probe assignment used to build the ARACNe-AP input (12,176 genes; see below) |
| `collectri_regulons_human.tsv` | CollecTRI TF–target table from OmniPath (`https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=1`, academic license; complexes removed, `mor = -1` for inhibition-only edges) |

The consolidated network is shipped with the repository so that the MRA can be run without Java or a large server; `run_aracne=true` regenerates it from `results/networks/SJS_matrix.txt` (about 10 min with `--cores 75`: each bootstrap is declared as 3 threads, since the JVM uses ~2.2 CPUs, and 4 GB RAM). The regenerated network has the same edges and p-values as the shipped one; MI values agree to ~1e-15 (order of summation across bootstraps).

**ARACNe-AP input and gene names.** The ARACNe-AP input matrix (190 pSS cases) comes from the original 2019 preprocessing of GSE51092, in which each Illumina probe was assigned to an HGNC symbol current at that time (likely via g:Profiler g:Convert); the exact assignment rule could not be recovered. `probe_map_GSE51092_2019.tsv` records that assignment, and `prep_GSE51092.R` rebuilds the input byte for byte from the GEO series matrix with it: every value is an unmodified GEO probe value. The expression matrix used for the signature, differential expression and MRA (`eset`) instead collapses probes by the platform `Symbol` annotation (11,409 genes). Because both use different gene names, `aracne2regulon()` keeps only network edges between genes named identically in both (7,900 of the 12,136 genes in the network).

### Steps

| Script | Purpose | Main outputs |
|--------|---------|--------------|
| `tf_list_lambert.R` | Lambert et al. 2018 TF list (1,639 TFs) | `results/networks/TFs_Lambert2018.txt` |
| `prep_GSE51092.R` | GSE51092 probe collapse and labels; cases-only ARACNe-AP input from the 2019 probe map | `results/mra/eset_GSE51092.rda`, `results/networks/SJS_matrix.txt` |
| `mra_GSE51092.R` | msVIPER with the ARACNe-AP regulon; limma | `results/mra/mra_GSE51092.rda`, `results/de/DE_limma_GSE51092.txt`, `results/figures/Top10_MRs_GSE51092.pdf` |
| `mra_GSE84844.R` | corto network, msVIPER; limma | `results/networks/regulon_GSE84844_corto.rda`, `results/mra/mra_GSE84844.rda`, `results/de/DE_limma_GSE84844.txt` |
| `compare_MRs.R` | Shared-TF table, RRA over the union of regulators | `results/integration/MRs_comparacion.txt`, `comparison_results.rda` |
| `meta_MRs.R` | Stouffer + restricted RRA over shared TFs (consensus MRs) | `results/integration/MRs_meta_shared.txt`, `MRs_concordance_tests.txt` |
| `figure_MRs.R` | Figure 3 | `results/figures/MRs_integrative_figure.pdf` |
| `ora_regulons.R` | ORA of consensus regulon targets (Figure 4) | `results/enrichment/ORA_*.txt`, `results/figures/ORA_combined.pdf` |
| `ora_nonconsensus.R` | Same ORA for HESX1, ZNF684, THRA, MEF2D and FOXJ2 (outside the consensus, evaluated with the fixed network) | `results/enrichment/ORA_nonconsensus.tsv` |
| `compare_ED.R` | Replicated DE genes, ORA | `results/de/DE_replicated_{up,down}.txt`, `comparison_DE_results.rda` |
| `figure_DE.R` | Figure 5 | `results/figures/DE_integrative_figure.pdf` |
| `meta_DE.R` | Stouffer + restricted RRA over shared genes | `results/integration/DE_meta_shared.txt` |
| `convergence_MR_DE.R` | Consensus MRs that are also replicated DE genes | `results/integration/MR_DE_convergence{,_summary}.txt` |
| `de_GSE66795.R` | Replication of DE genes in GSE66795 (Table S1) | `results/validation/Table_S1_DE_GSE66795_replication.xlsx` |
| `ssgsea_GSE140161.R` | ssGSEA of the replicated up signature in GSE140161 | `results/validation/ssGSEA_GSE140161_scores_replicated.txt`, `results/figures/ssGSEA_GSE140161_replicated.pdf` |
| `viper_cross_network.R` | Robustness: VIPER with the other cohort's network | `results/robustness/viper_cross_network*.txt` |
| `viper_collectri.R` | Robustness: VIPER with CollecTRI regulons | `results/robustness/viper_collectri*.txt`, `consensus_vs_collectri.txt` |
| `deconvolution_mra.R` | Robustness: cell composition (xCell, MCP-counter; BH over all 16 cell types of both methods within each cohort) and composition-adjusted MRA | `results/robustness/deconv_composition_tests.txt`, `mra_adjusted_*.txt` |
| `composition_reported.R` | Composition figures cited in the Supplementary Results | `results/robustness/composition_reported.tsv` |
| `fgsea_regulons.R` | Robustness: fgsea of regulons vs VIPER | `results/enrichment/fgsea_regulons_*.txt`, `fgsea_vs_viper_concordance.txt` |
| `table1_mr_summary.R` | Table 1: the 14 consensus MRs with activity, own-transcript DE, top ORA term and CollecTRI coverage | `results/integration/Table1_MR_summary.tsv` |
| `supplementary_tables.R` | Workbook of Supplementary Tables S1–S7 | `results/Supplementary_Tables_S1-S7.xlsx` |

`scripts/helpers.R` holds the shared probe-collapsing and MR-table functions. The full input/output graph is in the `Snakefile` (`snakemake --dag | dot -Tpdf > dag.pdf`).

---

## Output layout (`results/`, git-ignored)

| Folder | Contents |
|--------|----------|
| `results/networks/` | TF list, cases-only GSE51092 matrix (`SJS_matrix.txt`), ARACNE-AP network (`sjs_network.txt`; bootstraps in `aracne/` when re-inferred), corto regulon (`regulon_GSE84844_corto.rda`) |
| `results/mra/` | `eset_GSE51092.rda`, `mra_GSE51092.rda`, `mra_GSE84844.rda`, `top100_MRs_*.txt` |
| `results/de/` | Per-cohort limma tables, shared-gene comparison, RRA lists, replicated gene lists, `comparison_DE_results.rda` |
| `results/integration/` | Cross-cohort MR integration (`MRs_meta_shared.txt`, `MRs_concordance_tests.txt`, `comparison_results.rda`, `meta_MRs_results.rda`), DE integration (`DE_meta_shared.txt`, `meta_DE_results.rda`), MR/DE convergence tables |
| `results/enrichment/` | ORA tables and objects for MR regulons (`ORA_*`) and for replicated DE genes (`DE_ORA_*`) |
| `results/validation/` | GSE66795 replication (limma table, `.rda`, Supplementary Table S1 `.xlsx`) and GSE140161 ssGSEA scores |
| `results/robustness/` | Cross-network VIPER, CollecTRI VIPER, composition-adjusted MRA (steps 12–14) |
| `results/figures/` | All PDF figures |
| `results/Supplementary_Tables_S1-S7.xlsx` | Supplementary Tables S1–S7, one sheet per table plus an index |

### Key objects

| File | Contents |
|------|----------|
| `results/mra/mra_GSE51092.rda` | `mrs_gse51092`, `regulon_gse51092`, `eset`, `disease`, `de_gse51092`, `top_mrs` |
| `results/mra/mra_GSE84844.rda` | `mrs_gse84844`, `de_gse84844`, `top_mrs`, `eset`, `disease` |
| `results/integration/comparison_results.rda` | `shared` (TFs with a regulon in both networks), RRA results (`rra_act`, `rra_rep`) |
| `results/integration/meta_MRs_results.rda` | Stouffer + restricted RRA over shared TFs (`sh`, `concordance`) |
| `results/de/comparison_DE_results.rda` | `shared` genes with Stouffer z, replicated sets (`sig_up`, `sig_dn`), ORA objects |
| `results/integration/meta_DE_results.rda` | Stouffer + restricted RRA over shared genes (`d`) |
| `results/figures/MRs_integrative_figure.pdf` | Figure 3 — MR cross-dataset comparison |
| `results/figures/ORA_combined.pdf` | Figure 4 — ORA on consensus regulon target genes |
| `results/figures/DE_integrative_figure.pdf` | Figure 5 — replicated DE + ORA |

---

## Citation

Perez-Calixto D, Hernández-Lemus E, Tovar H. *Master Regulator Analysis Reveals a Reproducible Interferon-Centered Regulatory Architecture in Primary Sjögren's Syndrome.* Frontiers in Immunology (under review).

To cite this code, use the metadata in [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button) or the Zenodo DOI of the release.

---

## Data sources

- GSE51092: Becker et al. 2014 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092>
- GSE84844: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844>
- GSE66795: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66795>
- GSE140161: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140161>
- Lambert TF list: <https://humantfs.ccbr.utoronto.ca>
