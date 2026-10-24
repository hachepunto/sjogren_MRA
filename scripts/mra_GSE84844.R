library(GEOquery)
library(limma)
library(corto)
library(viper)
library(tibble)
library(readr)

source("scripts/helpers.R")

invisible(lapply(c("results/networks", "results/mra", "results/de"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

lambert_symbols <- readLines("results/networks/TFs_Lambert2018.txt")   # tf_list_lambert.R

# ── 1. Expression matrix (anotación GPL570 cacheada en data/) ────────────────
gse  <- getGEO(filename = "data/GSE84844_series_matrix.txt",
               GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data")
eset <- collapse_probes(exprs(gse), fData(gse)[["Gene symbol"]])
cat("Expression matrix:", nrow(eset), "genes x", ncol(eset), "samples\n")

# ── 2. Case / control labels ──────────────────────────────────────────────────
pdata   <- pData(gse)
disease <- ifelse(pdata$"disease:ch1" == "Healthy control", "control", "pSS")
print(table(disease))

# ── 3. TF hub list (Lambert, filtered to genes present in GSE84844) ───────────
centroids <- intersect(lambert_symbols, rownames(eset))
cat("TFs present in GSE84844:", length(centroids), "of", length(lambert_symbols), "\n")

# ── 4. Transcriptional network inference (corto) ──────────────────────────────
# Hilos como primer argumento (Snakemake pasa {threads}); cada bootstrap fija su
# propia semilla en el worker, así que la red no depende del número de hilos
N_THREADS <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(N_THREADS)) N_THREADS <- 4

set.seed(1)   # bootstraps de corto y permutaciones de ttestNull
regulon_gse84844 <- corto(
    inmat       = eset,
    centroids   = centroids,
    nbootstraps = 100,
    p           = 1e-8,
    nthreads    = N_THREADS,
    verbose     = TRUE
)

save(regulon_gse84844, disease, file = "results/networks/regulon_GSE84844_corto.rda")
cat("Network saved to results/networks/regulon_GSE84844_corto.rda\n")

# ── 5. Master Regulator Analysis — GSE84844 ───────────────────────────────────
Ncases    <- which(disease == "pSS")
Ncontrols <- which(disease == "control")

signature <- rowTtest(eset[, Ncases], eset[, Ncontrols])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]

nullmodel <- ttestNull(eset[, Ncases], eset[, Ncontrols], per = 1000, repos = TRUE)

mrs_gse84844 <- msviper(signature, regulon_gse84844, nullmodel)

print(summary(mrs_gse84844))

top_mrs <- viper_mrsTopTable(mrs_gse84844, n = 100)
write_tsv(top_mrs, "results/mra/top100_MRs_GSE84844.txt")

# ── 6. Expresión diferencial (limma) ──────────────────────────────────────────
disease_f    <- factor(disease, levels = c("control", "pSS"))
design       <- model.matrix(~ disease_f)
fit          <- lmFit(eset, design)
fit          <- eBayes(fit)
de_gse84844 <- topTable(fit, coef = 2, number = Inf,
                         adjust.method = "BH", sort.by = "P")
de_gse84844 <- tibble::rownames_to_column(de_gse84844, var = "gene")
write_tsv(de_gse84844, "results/de/DE_limma_GSE84844.txt")
cat("DE: ", sum(de_gse84844$adj.P.Val < 0.05), "genes significativos (FDR < 0.05)\n")

save(mrs_gse84844, top_mrs, de_gse84844, eset, disease, file = "results/mra/mra_GSE84844.rda")

