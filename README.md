# Master Regulator Analysis of Primary Sjögren's Syndrome

Code for the article *Master Regulator Analysis Reveals a Reproducible Interferon-Centered Regulatory Architecture in Primary Sjögren's Syndrome* (Perez-Calixto, Hernández-Lemus and Tovar; Frontiers in Immunology, under review).

The pipeline infers a transcriptional regulatory network in each of two independent whole-blood cohorts of primary Sjögren's syndrome (pSS), estimates transcription factor (TF) activity with VIPER, and integrates both cohorts over the regulators and genes measured in both. It then characterizes the consensus regulators functionally, validates the replicated expression signature in two further cohorts, and tests the robustness of the consensus to network inference, interactome source and cell composition. The whole analysis is a [Snakemake](https://snakemake.readthedocs.io) workflow that goes from the public GEO files to every figure and table of the article.

---

## Analysis overview

1. **Regulatory networks.** ARACNe-AP on the 190 pSS cases of GSE51092 and corto on GSE84844, both with the 1,639 human TFs of Lambert et al. 2018 as candidate regulators.
2. **Master regulator analysis.** msVIPER on the pSS-versus-control signature of each cohort, with its own network and a null model of 1,000 sample permutations.
3. **Integration of regulator activity.** The NES of the two cohorts are combined by Stouffer's method, z = (NES₁ + NES₂)/√2, over the 240 TFs with a regulon in both networks, and corrected by Benjamini-Hochberg. **Consensus master regulators: same sign in both cohorts and BH FDR < 0.10** (14 regulators: 12 activated, 2 repressed).
4. **Differential expression.** limma in each cohort. **Replicated genes: FDR < 0.05 in each cohort with the same sign** (375 up, 200 down, over 8,007 shared genes). A Stouffer meta-analysis of the limma p-values is reported as a complementary summary.
5. **Functional characterization.** Over-representation analysis (GO Biological Process and Reactome) of the consensus regulons and of the replicated genes; Table 1 summarizes each consensus regulator.
6. **External validation.** The replicated genes are tested in GSE66795 (161 of the 375 up- and 68 of the 200 down-regulated genes recover with the same sign) and the replicated up-regulated signature is scored by ssGSEA in the 351 pSS patients of GSE140161.
7. **Robustness.** VIPER with the network of the other cohort, VIPER with curated CollecTRI regulons, master regulator analysis adjusted for cell composition (xCell, MCP-counter), and rank-based enrichment of the regulons with fgsea.

Robust Rank Aggregation (RRA), the integration criterion of the originally submitted manuscript (v1.x of this repository), is still computed by `compare_MRs.R`, `compare_ED.R`, `meta_MRs.R` and `meta_DE.R`, but only as a descriptive reference; it does not define any result.

---

## Cohorts

| | Discovery | Replication |
|---|---|---|
| **GEO accession** | [GSE51092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092) | [GSE84844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844) |
| **Platform** | Illumina HumanWG-6 v3.0 (GPL6884) | Affymetrix Human Genome U133 Plus 2.0 (GPL570) |
| **Design** | 32 controls + 190 pSS (222) | 30 controls + 30 pSS (60) |
| **Network** | ARACNe-AP (190 pSS cases) | corto (60 samples) |

External validation cohorts (not used for network inference or for the consensus):

- [GSE66795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66795): 29 controls + 131 pSS, Illumina HumanHT-12 v4.0 (GPL10558); validation of the replicated genes (Table S1).
- [GSE140161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140161): 351 pSS, no controls, Affymetrix Clariom S (GPL23159); ssGSEA of the replicated up-regulated signature.

All four cohorts are used as the preprocessed, submitter-normalized data deposited in GEO.

---

## Quick start

```bash
git clone --branch v2.0.0 https://github.com/hachepunto/sjogren_MRA.git   # version used for the article
cd sjogren_MRA
# R 4.6 environment with the packages of envs/install_r_packages.R (see Setup)
snakemake --cores 8
```

With the shipped ARACNe-AP network the full workflow takes about 20 minutes on 8 cores once the inputs are downloaded (about 300 MB in `data/`). Results go to `results/` and one log per step to `logs/`.

---

## Setup

### R environment

R 4.6.1 with Bioconductor 3.23. `envs/install_r_packages.R` installs every package the pipeline needs: CRAN and Bioconductor packages, corto 1.2.4, and xCell and MCPcounter from GitHub pinned to the commits used for the article. `package_versions.tsv` lists the exact versions used for the article, including the annotation databases (org.Hs.eg.db, GO.db, reactome.db) that determine the enrichment results; `envs/package_versions.R` regenerates it from the active environment. With conda, for example:

```bash
conda create -n sjogren -c conda-forge r-base=4.6.1 c-compiler cxx-compiler fortran-compiler make cmake pkg-config \
    libcurl libxml2 openssl zlib xz bzip2 icu glpk gmp freetype fontconfig harfbuzz fribidi libpng libjpeg-turbo libtiff cairo libuv imagemagick
conda activate sjogren      # packages must compile with the environment's own toolchain
Rscript envs/install_r_packages.R
```

Bioconductor packages are fixed by the Bioconductor release; CRAN packages are installed at their current versions, so they may differ from `package_versions.tsv` unless a dated CRAN snapshot is used (e.g. `options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-10-20"))` before running the script).

The workflow itself needs Snakemake (≥ 8; tested with 9.23.1), `curl`, `gunzip` and `md5sum`. Re-inferring the network also needs Java and the [ARACNe-AP](https://github.com/califano-lab/ARACNe-AP) jar (built with `ant main`); the re-inference was tested with a jar built in August 2023 (MD5 `1292fe68e1336649f4d17410e22e8f2e`) run with OpenJDK 17.

### Configuration (`config.yaml`)

| Key | Default | Meaning |
|-----|---------|---------|
| `rscript` | `Rscript` | Path to the `Rscript` of the R 4.6 installation |
| `corto_threads` | `4` | Threads for corto (GSE84844 network) |
| `run_aracne` | `false` | `false`: use the shipped network in `resources/`; `true`: re-infer it with ARACNe-AP |
| `aracne_jar` | placeholder path | Path to `aracne.jar` (only with `run_aracne: true`) |
| `aracne_bootstraps`, `aracne_pvalue`, `aracne_mem` | `100`, `1E-8`, `4G` | ARACNe-AP parameters; each bootstrap is one job with `aracne_mem` of Java heap |

Any key can be overridden on the command line, e.g. `--config rscript=/opt/R/4.6.1/bin/Rscript`.

---

## Running

```bash
snakemake --cores 8                                                             # full workflow, shipped network
snakemake --cores 74 --config run_aracne=true aracne_jar=/path/to/aracne.jar    # re-infer the ARACNe-AP network
snakemake -n                                                                    # dry run: list pending jobs
```

- **Downloads.** Snakemake downloads the GEO series matrices, the GPL platform annotations that GEOquery uses to map probes to genes, the GSE140161 processed matrix, the Lambert et al. TF table and the MCP-counter gene signature (from the pinned MCPcounter commit) into `data/`. `results/data_checksums.txt` compares every download with the MD5 sums of the files used for the article (`resources/data_md5.txt`) and flags any file that GEO has updated since.
- **CPU use.** Every rule declares the number of CPUs it uses in `threads:`, and BLAS, OpenMP, readr and data.table are limited to one thread per process, so `--cores` bounds the CPU use of the whole workflow. On a shared server it can also be enforced with CPU affinity, e.g. `taskset -c 0-73 snakemake --cores 74`.
- **Re-running.** Each rule lists its script among its inputs, so editing a script re-runs that step and everything downstream. Every script can also be run on its own from the repository root (`Rscript scripts/<name>.R`) once its inputs exist.

### Re-inferring the ARACNe-AP network

`run_aracne=true` computes the mutual-information threshold, runs 100 bootstraps (seeds 1–100) as independent jobs and consolidates them. Each bootstrap is declared as 3 threads (a JVM uses about 2.2 CPUs) and 5 GB of RAM; Snakemake enforces the memory only if a limit is given, e.g. `--resources mem_mb=120000`. With `--cores 74` the network takes about 6 minutes and the whole workflow, downloads included, about 27 minutes. The re-inferred network has the same 140,121 edges and p-values as the shipped one, and mutual information agrees to within 1e-15 (the differences come from the order of summation across bootstraps).

---

## Shipped resources (`resources/`)

| File | Contents |
|------|----------|
| `sjs_network_aracne.txt.gz` | Consolidated ARACNe-AP network of GSE51092 used in the article (190 pSS cases, Lambert TFs, p < 1e-8, 100 bootstraps with seeds 1–100) |
| `probe_map_GSE51092_2019.tsv` | Gene → Illumina probe assignment used to build the ARACNe-AP input (12,176 genes; see below) |
| `collectri_regulons_human.tsv` | CollecTRI TF–target interactions (41,212 interactions, 1,176 TFs) downloaded from the [OmniPath](https://omnipathdb.org) web service on 22 October 2026; complexes removed, `mor = -1` for interactions annotated only as inhibitory and `+1` otherwise. `scripts/get_collectri.R` rebuilds it (outside the workflow); on 24 October 2026 it reproduced the shipped file byte for byte |
| `data_md5.txt` | MD5 sums of the downloaded input files used for the article |

**ARACNe-AP input and gene names.** The ARACNe-AP input matrix comes from the original 2019 preprocessing of GSE51092, in which Illumina probes were mapped to the HGNC symbols current at that time (most likely with g:Profiler g:Convert) and one probe was kept per gene; the exact assignment rule could not be recovered. `probe_map_GSE51092_2019.tsv` records that assignment, and `prep_GSE51092.R` rebuilds the input byte for byte from the GEO series matrix with it, so every value is an unmodified GEO probe value. The expression matrix used for the signature, differential expression and master regulator analysis instead collapses probes by the platform `Symbol` annotation (11,409 genes). Because the two use different gene names, `aracne2regulon()` keeps only the network edges between genes named identically in both (7,900 of the 12,136 genes in the network).

**CollecTRI license.** CollecTRI is assembled by OmniPath from several resources, each with its own license; the file is redistributed here only so that the analysis can be reproduced, under the terms of those original sources (see [OmniPath licenses](https://omnipathdb.org/info)). It is not covered by the GPL-3.0 license of this repository.

---

## Workflow steps

| Script | Purpose | Main outputs |
|--------|---------|--------------|
| **Networks and per-cohort analysis** | | |
| `tf_list_lambert.R` | Lambert et al. 2018 TF list (1,639 TFs) | `results/networks/TFs_Lambert2018.txt` |
| `prep_GSE51092.R` | GSE51092 probe collapse and labels; cases-only ARACNe-AP input from the 2019 probe map | `results/mra/eset_GSE51092.rda`, `results/networks/SJS_matrix.txt` |
| `mra_GSE51092.R` | msVIPER with the ARACNe-AP regulon; limma | `results/mra/mra_GSE51092.rda`, `results/de/DE_limma_GSE51092.txt`, `results/figures/Top10_MRs_GSE51092.pdf` |
| `mra_GSE84844.R` | corto network, msVIPER; limma | `results/networks/regulon_GSE84844_corto.rda`, `results/mra/mra_GSE84844.rda`, `results/de/DE_limma_GSE84844.txt` |
| **Integration of regulator activity** | | |
| `compare_MRs.R` | Table of the 240 shared TFs; NES scatter; RRA over the union of regulators (reference only) | `results/integration/MRs_comparacion.txt`, `MRs_RRA_*.txt`, `comparison_results.rda`, `results/figures/MRs_scatter_NES.pdf` |
| `meta_MRs.R` | Stouffer meta-analysis over the shared TFs; **consensus master regulators**; restricted RRA (descriptive) (Table S2) | `results/integration/MRs_meta_shared.txt` (column `consensus`), `MRs_concordance_tests.txt` |
| `figure_MRs.R` | Figure 3 | `results/figures/MRs_integrative_figure.pdf` |
| `ora_regulons.R` | ORA of the ARACNe-AP regulons of the 14 consensus regulators (Figure 4) | `results/enrichment/ORA_{GO_BP,Reactome}.txt`, `results/figures/ORA_combined.pdf` |
| `ora_nonconsensus.R` | Same ORA for the five regulators outside the consensus that the article discusses: HESX1, ZNF684 and THRA (no corto regulon of at least 25 targets) and MEF2D and FOXJ2 (concordantly repressed, just above the FDR threshold) | `results/enrichment/ORA_nonconsensus.tsv` |
| **Differential expression** | | |
| `compare_ED.R` | Shared genes; **replicated genes**; ORA of the replicated genes; RRA over the full tables (reference only) | `results/de/DE_comparacion.txt`, `DE_replicated_{up,down}.txt`, `comparison_DE_results.rda`, `results/enrichment/DE_ORA_*.txt` |
| `figure_DE.R` | Figure 5 | `results/figures/DE_integrative_figure.pdf` |
| `meta_DE.R` | Stouffer meta-analysis over the 8,007 shared genes; restricted RRA (descriptive) (Table S3) | `results/integration/DE_meta_shared.txt` |
| `convergence_MR_DE.R` | Consensus regulators whose own transcript is a replicated gene | `results/integration/MR_DE_convergence{,_summary}.txt` |
| `table1_mr_summary.R` | Table 1: activity, own-transcript expression, top enriched term and CollecTRI coverage of each consensus regulator | `results/integration/Table1_MR_summary.tsv` |
| **Validation** | | |
| `de_GSE66795.R` | Validation of the replicated genes in GSE66795 (Table S1) | `results/validation/Table_S1_DE_GSE66795_replication.xlsx` |
| `ssgsea_GSE140161.R` | ssGSEA of the replicated up-regulated signature in GSE140161 (Figure S1) | `results/validation/ssGSEA_GSE140161_scores_replicated.txt`, `results/figures/ssGSEA_GSE140161_replicated.pdf` |
| **Robustness** | | |
| `viper_cross_network.R` | VIPER with the network of the other cohort and fixed-network meta-analysis (Table S4) | `results/robustness/viper_cross_network*.txt` |
| `viper_collectri.R` | VIPER with CollecTRI regulons (Table S5) | `results/robustness/viper_collectri*.txt`, `consensus_vs_collectri.txt` |
| `deconvolution_mra.R` | Cell composition with xCell and MCP-counter (BH over the 16 cell types of both methods within each cohort) and MRA adjusted for six xCell populations (Table S6) | `results/robustness/deconv_composition_tests.txt`, `{xcell,mcpcounter}_*.txt`, `mra_adjusted_*.txt` |
| `composition_reported.R` | Composition results cited in the Supplementary Results | `results/robustness/composition_reported.tsv` |
| `fgsea_regulons.R` | fgsea of the positive- and negative-mode targets of every regulon versus VIPER (Table S7) | `results/enrichment/fgsea_regulons_*.txt`, `fgsea_vs_viper_concordance.txt` |
| **Supplementary tables** | | |
| `supplementary_tables.R` | Workbook of Supplementary Tables S1–S7 | `results/Supplementary_Tables_S1-S7.xlsx` |

`scripts/helpers.R` holds the functions shared across scripts (probe collapsing, signature, msVIPER call, symbol mapping, consensus list, MR tables). The full input/output graph is in the `Snakefile` (`snakemake --dag | dot -Tpdf > dag.pdf`).

---

## Article figures and tables

| Item | File |
|------|------|
| Figure 1 (workflow schematic) | Drawn in BioRender, not generated by the code |
| Figure 2 | `results/figures/Top10_MRs_GSE51092.pdf` |
| Figure 3 | `results/figures/MRs_integrative_figure.pdf` |
| Figure 4 | `results/figures/ORA_combined.pdf` |
| Figure 5 | `results/figures/DE_integrative_figure.pdf` |
| Figure 6 (conceptual model) | Conceptual diagram, not generated by the code |
| Table 1 | `results/integration/Table1_MR_summary.tsv` |
| Supplementary Tables S1–S7 | `results/Supplementary_Tables_S1-S7.xlsx` |
| Supplementary Figure S1 | `results/figures/ssGSEA_GSE140161_replicated.pdf` |

---

## Output layout (`results/`, git-ignored)

| Folder | Contents |
|--------|----------|
| `results/networks/` | TF list, cases-only GSE51092 matrix (`SJS_matrix.txt`), ARACNe-AP network (`sjs_network.txt`; threshold and bootstraps in `aracne/` when re-inferred), corto regulon (`regulon_GSE84844_corto.rda`) |
| `results/mra/` | `eset_GSE51092.rda`, `mra_GSE51092.rda`, `mra_GSE84844.rda`, `top100_MRs_*.txt` (significant regulators, p < 0.05, up to 100) |
| `results/de/` | Per-cohort limma tables, shared-gene table, replicated gene lists, RRA reference lists, `comparison_DE_results.rda` |
| `results/integration/` | Regulator integration (`MRs_comparacion.txt`, `MRs_meta_shared.txt`, `MRs_concordance_tests.txt`, `MRs_RRA_*.txt`, `comparison_results.rda`, `meta_MRs_results.rda`), expression integration (`DE_meta_shared.txt`, `meta_DE_results.rda`), MR/DE convergence and `Table1_MR_summary.tsv` |
| `results/enrichment/` | ORA of the consensus regulons (`ORA_*`), of the five non-consensus regulons (`ORA_nonconsensus.*`) and of the replicated genes (`DE_ORA_*`); fgsea of the regulons (`fgsea_*`) |
| `results/validation/` | GSE66795 replication (limma table, `.rda`, Table S1 `.xlsx`) and GSE140161 ssGSEA scores |
| `results/robustness/` | Cross-network VIPER, CollecTRI VIPER, deconvolution scores and tests, composition-adjusted MRA |
| `results/figures/` | All PDF figures |
| `results/Supplementary_Tables_S1-S7.xlsx` | Supplementary Tables S1–S7, one sheet per table plus an index |
| `results/data_checksums.txt` | MD5 check of the downloaded inputs |

### Key objects

| File | Contents |
|------|----------|
| `results/mra/eset_GSE51092.rda` | `eset` (11,409 genes × 222 samples), `disease` |
| `results/mra/mra_GSE51092.rda` | `mrs_gse51092`, `regulon_gse51092`, `eset`, `disease`, `de_gse51092`, `top_mrs` |
| `results/networks/regulon_GSE84844_corto.rda` | `regulon_gse84844` (corto regulon), `disease` |
| `results/mra/mra_GSE84844.rda` | `mrs_gse84844`, `de_gse84844`, `top_mrs`, `eset`, `disease` |
| `results/integration/comparison_results.rda` | `sum_gse51092`, `sum_84844` (all regulators per cohort), `shared` (the 240 shared TFs), `rra_act`, `rra_rep` (RRA over the union, reference) |
| `results/integration/meta_MRs_results.rda` | `sh` (Stouffer meta-analysis and consensus), `concordance`, `rra_shared_act`, `rra_shared_rep` (restricted RRA) |
| `results/de/comparison_DE_results.rda` | `shared` genes with Stouffer z, replicated sets (`sig_up`, `sig_dn`), RRA over the full tables (`rra_up`, `rra_dn`), ORA objects (`ego_de`, `ereact_de`) |
| `results/integration/meta_DE_results.rda` | `d` (Stouffer meta-analysis and replication over shared genes), `rra_shared_up`, `rra_shared_dn` |

---

## Reproducibility notes

- **Determinism.** The permutation null models are seeded (`set.seed(1)`) and corto seeds each of its 100 bootstraps internally, so its network does not depend on the number of threads. A complete re-run from the GEO downloads, including the re-inference of the ARACNe-AP network, reproduces the network and every table of the article up to floating-point noise (≤ 1e-14 in NES and p-values; the order of genes inside the `geneID` strings of the ORA tables may also change).
- **Gene identifiers.** Enrichment analyses map gene symbols to Entrez IDs with org.Hs.eg.db. About 30% of the GSE51092 platform symbols do not map (mostly outdated aliases and LOC/FLJ/DKFZ placeholders); they are removed from both the gene lists and the background, so the tests are conditional on the mappable genes and lose some power.
- **Ties in enrichment results.** Terms with identical statistics (for example Reactome pathways that share the same proteasome subunits in the ZNF227 regulon) can be listed in either order, which can change which term is drawn last in a panel of Figure 4.

---

## Citation

Perez-Calixto D, Hernández-Lemus E, Tovar H. *Master Regulator Analysis Reveals a Reproducible Interferon-Centered Regulatory Architecture in Primary Sjögren's Syndrome.* Frontiers in Immunology (under review).

To cite this code, use the metadata in [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button) or the Zenodo DOI of the release.

## License

The code is released under the [GNU General Public License v3.0](LICENSE). The CollecTRI table in `resources/` keeps the licenses of its original sources (see above).

## Methods references

- ARACNe-AP: Lachmann A et al. *Bioinformatics* 2016;32:2233–2235. doi:[10.1093/bioinformatics/btw216](https://doi.org/10.1093/bioinformatics/btw216)
- corto: Mercatelli D et al. *Bioinformatics* 2020;36:3916–3917. doi:[10.1093/bioinformatics/btaa223](https://doi.org/10.1093/bioinformatics/btaa223)
- VIPER: Alvarez MJ et al. *Nat Genet* 2016;48:838–847. doi:[10.1038/ng.3593](https://doi.org/10.1038/ng.3593)
- limma: Ritchie ME et al. *Nucleic Acids Res* 2015;43:e47. doi:[10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)
- Robust Rank Aggregation: Kolde R et al. *Bioinformatics* 2012;28:573–580. doi:[10.1093/bioinformatics/btr709](https://doi.org/10.1093/bioinformatics/btr709)
- GSVA/ssGSEA: Hänzelmann S et al. *BMC Bioinformatics* 2013;14:7. doi:[10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)
- xCell: Aran D et al. *Genome Biol* 2017;18:220. doi:[10.1186/s13059-017-1349-1](https://doi.org/10.1186/s13059-017-1349-1)
- MCP-counter: Becht E et al. *Genome Biol* 2016;17:218. doi:[10.1186/s13059-016-1070-5](https://doi.org/10.1186/s13059-016-1070-5)
- fgsea: Korotkevich G et al. *bioRxiv* 2021. doi:[10.1101/060012](https://doi.org/10.1101/060012)

## Data sources

- GSE51092: Lessard CJ et al. Variants at multiple loci implicated in both innate and adaptive immune responses are associated with Sjögren's syndrome. *Nat Genet* 2013;45:1284–1292. doi:[10.1038/ng.2792](https://doi.org/10.1038/ng.2792)
- GSE84844: Tasaki S et al. Multiomic disease signatures converge to cytotoxic CD8 T cells in primary Sjögren's syndrome. *Ann Rheum Dis* 2017;76:1458–1466. doi:[10.1136/annrheumdis-2016-210788](https://doi.org/10.1136/annrheumdis-2016-210788)
- GSE66795: James K et al. A transcriptional signature of fatigue derived from patients with primary Sjögren's syndrome. *PLoS One* 2015;10:e0143970. doi:[10.1371/journal.pone.0143970](https://doi.org/10.1371/journal.pone.0143970)
- GSE140161: Felten R et al. Identification of new candidate drugs for primary Sjögren's syndrome using a drug repurposing transcriptomic approach. *Rheumatology* 2023;62:3715–3723. doi:[10.1093/rheumatology/kead096](https://doi.org/10.1093/rheumatology/kead096)
- Lambert SA et al. The human transcription factors. *Cell* 2018;172:650–665. doi:[10.1016/j.cell.2018.01.029](https://doi.org/10.1016/j.cell.2018.01.029); <https://humantfs.ccbr.utoronto.ca>.
- CollecTRI: Müller-Dott S et al. *Nucleic Acids Res* 2023;51:10934–10949. doi:[10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841); distributed through OmniPath (Türei D et al. *Mol Syst Biol* 2021;17:e9923. doi:[10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923)).
