# Instala los paquetes de R del pipeline (R 4.6 / Bioconductor 3.23).
# Correr con el ambiente de R activado para que compile con sus propios compiladores.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)
BiocManager::install(c(
    "readr", "tibble", "purrr", "ggplot2", "patchwork", "openxlsx", "remotes",
    "RobustRankAggreg", "quadprog",
    "GEOquery", "matrixStats", "limma", "viper", "fgsea",
    "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "GSVA"
), ask = FALSE, update = FALSE)
# corto está en CRAN; se fija la versión usada en el manuscrito
remotes::install_version("corto", version = "1.2.4", upgrade = "never")
remotes::install_github(c("dviraran/xCell", "ebecht/MCPcounter/Source"),
                        upgrade = "never")
