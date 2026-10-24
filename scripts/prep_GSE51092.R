library(GEOquery)

source("scripts/helpers.R")

invisible(lapply(c("results/networks", "results/mra"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ── 1. Expression matrix completa (15063 probes → 11409 genes x 222 muestras) ─
gse  <- getGEO(filename = "data/GSE51092_series_matrix.txt", GSEMatrix = TRUE, destdir = "data")
eset <- collapse_probes(exprs(gse), fData(gse)[["Symbol"]])
cat("Expression matrix:", nrow(eset), "genes x", ncol(eset), "samples\n")

# ── 2. Disease labels ─────────────────────────────────────────────────────────
series    <- readLines("data/GSE51092_series_matrix.txt")
char_line <- grep("disease", series, value = TRUE)[1]
labels    <- gsub('"', '', unlist(strsplit(char_line, "\t"))[-1])
disease   <- ifelse(labels == "disease state: none", "control", "pSS")
print(table(disease))

save(eset, disease, file = "results/mra/eset_GSE51092.rda")

# ── 3. Matriz solo de casos para ARACNE-AP ───────────────────────────────────
# Usa la asignación gen → sonda del preprocesamiento original (2019; símbolos
# HGNC de esa fecha), distinta del colapso por "Symbol" del paso 1. Reproduce
# byte a byte la entrada con la que se infirió la red del manuscrito.
probe_map <- read.delim("resources/probe_map_GSE51092_2019.tsv")
cases_mat <- exprs(gse)[probe_map$probe, disease == "pSS"]
colnames(cases_mat) <- paste0("case_", seq_len(ncol(cases_mat)))
write.table(data.frame(genes = probe_map$gene, cases_mat, check.names = FALSE),
            "results/networks/SJS_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Matriz ARACNE-AP:", nrow(cases_mat), "genes x", ncol(cases_mat), "casos\n")
