library(GEOquery)
library(matrixStats)
library(limma)
library(corto)
library(viper)
library(tibble)
library(purrr)
library(readr)

viper_mrsTopTable <- function(mrs, n = sum(mrs$es$p.value < p_threshold), p_threshold = 0.05) {
    sig <- mrs$es$p.value < p_threshold
    mrs_tops <- tibble(
        TF      = names(mrs$es$nes[sig]),
        nes     = mrs$es$nes[sig],
        p.value = mrs$es$p.value[sig]
    )
    mrs_tops <- mrs_tops[order(abs(mrs_tops$nes), decreasing = TRUE), ]
    mrs_tops <- mrs_tops[1:n, ]
    if (!is.null(mrs$ledge)) {
        mrs_tops$leading_edge <- sapply(mrs_tops$TF, function(tf) {
            le <- mrs$ledge[[tf]]
            if (is.null(le)) NA_character_ else paste(le, collapse = ",")
        })
    }
    return(mrs_tops)
}

packinfo <- function(save = FALSE) {
    lst <- sessionInfo()$otherPkgs %>%
        purrr::map_df(~ data.frame(Package = .x$Package, Version = .x$Version))
    if (save) {
        write_delim(lst, "package_versions.tsv", delim = "\t")
    } else {
        return(lst)
    }
}

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# ── 0. TF list — Lambert et al. 2018 (Cell, PMID: 29425488) ──────────────────
# Source: https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
lambert_csv <- "data/lambert2018_DatabaseExtract_v1.01.csv"
if (!file.exists(lambert_csv)) {
    download.file(
        url      = "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv",
        destfile = lambert_csv,
        mode     = "wb"
    )
}
lambert_db      <- read.csv(lambert_csv, check.names = FALSE)
lambert_symbols <- unique(na.omit(lambert_db[["HGNC symbol"]][lambert_db[["Is TF?"]] == "Yes"]))
cat("Lambert et al. TFs:", length(lambert_symbols), "\n")

writeLines(lambert_symbols, "data/TFs_Lambert2018.txt")

# ── 1. Read local series matrix (GPL570 annotation downloaded once, cached) ──
gse <- getGEO(filename = "data/GSE84844_series_matrix.txt",
              GSEMatrix = TRUE,
              AnnotGPL  = TRUE)

# ── 2. Expression matrix & gene symbol mapping ────────────────────────────────
eset         <- exprs(gse)
fdata        <- fData(gse)
gene_symbols <- fdata[["Gene symbol"]]

keep         <- !is.na(gene_symbols) & gene_symbols != "" & !grepl("///", gene_symbols)
eset         <- eset[keep, ]
gene_symbols <- gene_symbols[keep]

# Collapse multiple probes per gene: keep the probe with highest median expression
eset <- do.call(rbind, lapply(
    split(seq_len(nrow(eset)), gene_symbols),
    function(idx) {
        if (length(idx) == 1) return(eset[idx, , drop = FALSE])
        eset[idx[which.max(rowMedians(eset[idx, , drop = FALSE]))], , drop = FALSE]
    }
))
rownames(eset) <- sort(unique(gene_symbols))
cat("Expression matrix:", nrow(eset), "genes x", ncol(eset), "samples\n")

# ── 3. Case / control labels ──────────────────────────────────────────────────
pdata   <- pData(gse)
disease <- ifelse(pdata$"disease:ch1" == "Healthy control", "control", "pSS")
print(table(disease))

# ── 4. TF hub list (Lambert, filtered to genes present in GSE84844) ───────────
centroids <- intersect(lambert_symbols, rownames(eset))
cat("TFs present in GSE84844:", length(centroids), "of", length(lambert_symbols), "\n")

# ── 5. Transcriptional network inference (corto) ──────────────────────────────
N_THREADS <- 4   # local run; set to 75 if running on server

regulon_gse84844 <- corto(
    inmat       = eset,
    centroids   = centroids,
    nbootstraps = 100,
    p           = 1e-8,
    nthreads    = N_THREADS,
    verbose     = TRUE
)

save(regulon_gse84844, disease, file = "regulon_GSE84844_corto.rda")
cat("Network saved to regulon_GSE84844_corto.rda\n")

# ── 6. Master Regulator Analysis — GSE84844 ───────────────────────────────────
Ncases    <- which(disease == "pSS")
Ncontrols <- which(disease == "control")

signature <- rowTtest(eset[, Ncases], eset[, Ncontrols])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]

nullmodel <- ttestNull(eset[, Ncases], eset[, Ncontrols], per = 1000, repos = TRUE)

mrs_gse84844 <- msviper(signature, regulon_gse84844, nullmodel)

print(summary(mrs_gse84844))

top_mrs <- viper_mrsTopTable(mrs_gse84844, n = 100)
write_tsv(top_mrs, "top100_MRs_GSE84844.txt")

# ── 7. Expresión diferencial (limma) ──────────────────────────────────────────
disease_f    <- factor(disease, levels = c("control", "pSS"))
design       <- model.matrix(~ disease_f)
fit          <- lmFit(eset, design)
fit          <- eBayes(fit)
de_gse84844 <- topTable(fit, coef = 2, number = Inf,
                         adjust.method = "BH", sort.by = "P")
de_gse84844 <- tibble::rownames_to_column(de_gse84844, var = "gene")
write_tsv(de_gse84844, "DE_limma_GSE84844.txt")
cat("DE: ", sum(de_gse84844$adj.P.Val < 0.05), "genes significativos (FDR < 0.05)\n")

save(mrs_gse84844, top_mrs, de_gse84844, file = "mra_GSE84844.rda")

packinfo(save = TRUE)
