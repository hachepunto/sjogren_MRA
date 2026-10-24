# Lista de TFs — Lambert et al. 2018 (Cell, PMID: 29425488)
# Fuente: https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)

lambert_db      <- read.csv("data/lambert2018_DatabaseExtract_v1.01.csv", check.names = FALSE)
lambert_symbols <- unique(na.omit(lambert_db[["HGNC symbol"]][lambert_db[["Is TF?"]] == "Yes"]))
cat("Lambert et al. TFs:", length(lambert_symbols), "\n")

writeLines(lambert_symbols, "results/networks/TFs_Lambert2018.txt")
