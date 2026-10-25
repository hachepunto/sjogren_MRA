# Instala los paquetes de R del pipeline (R 4.6 / Bioconductor 3.23).
# Correr con el ambiente de R activado para que compile con sus propios compiladores.
# Las versiones usadas para el artículo están en package_versions.tsv.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)
BiocManager::install(c(
    "readr", "tibble", "ggplot2", "ggrepel", "patchwork", "openxlsx", "remotes",
    "RobustRankAggreg", "quadprog",
    "GEOquery", "matrixStats", "limma", "viper", "fgsea",
    "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "GSVA"
), ask = FALSE, update = FALSE)
# corto está en CRAN; se fija la versión usada en el artículo
remotes::install_version("corto", version = "1.2.4", upgrade = "never")
# xCell y MCPcounter no tienen versiones publicadas: se fijan los commits usados
remotes::install_github(c("dviraran/xCell@20e2919eefd37e15af35f29f4944e30697098a28",
                          "ebecht/MCPcounter/Source@b6eac73e91c246fcff0bb1a5c68a816cd588fc48"),
                        upgrade = "never")
