library(GEOquery)
library(matrixStats)
library(limma)
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

# ── 1. Expression matrix completa (15063 probes → 11409 genes x 222 muestras) ─
gse   <- getGEO(filename = "data/GSE51092_series_matrix.txt", GSEMatrix = TRUE)
eset  <- exprs(gse)
symbols <- fData(gse)[["Symbol"]]

keep    <- !is.na(symbols) & symbols != "" & !grepl("///", symbols)
eset    <- eset[keep, ]
symbols <- symbols[keep]

eset <- do.call(rbind, lapply(
    split(seq_len(nrow(eset)), symbols),
    function(idx) {
        if (length(idx) == 1) return(eset[idx, , drop = FALSE])
        eset[idx[which.max(rowMedians(eset[idx, , drop = FALSE]))], , drop = FALSE]
    }
))
rownames(eset) <- sort(unique(symbols))
cat("Expression matrix:", nrow(eset), "genes x", ncol(eset), "samples\n")

# ── 2. Disease labels ─────────────────────────────────────────────────────────
series    <- readLines("data/GSE51092_series_matrix.txt")
char_line <- grep("disease", series, value = TRUE)[1]
labels    <- gsub('"', '', unlist(strsplit(char_line, "\t"))[-1])
disease   <- ifelse(labels == "disease state: none", "control", "pSS")
print(table(disease))

# ── 3. Export cases-only matrix for ARACNE-AP ────────────────────────────────
# Run run_aracne_ap.sh after this block, then re-open to continue from step 4.
# Skip if data/sjs_network.txt already exists.
cases_mat <- eset[, disease == "pSS"]
write.table(cases_mat, "data/SJS_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)

# ── 4. Regulon desde ARACNE-AP ────────────────────────────────────────────────
regulon_gse51092 <- aracne2regulon("data/sjs_network.txt", eset)
cat("Regulon:", length(regulon_gse51092), "TFs\n")

# ── 5. Firma molecular ────────────────────────────────────────────────────────
Ncases    <- which(disease == "pSS")
Ncontrols <- which(disease == "control")

signature <- rowTtest(eset[, Ncases], eset[, Ncontrols])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) *
              sign(signature$statistic))[, 1]

# ── 5. Null model ─────────────────────────────────────────────────────────────
nullmodel <- ttestNull(eset[, Ncases], eset[, Ncontrols],
                       per = 1000, repos = TRUE)

# ── 6. Master Regulator Analysis ──────────────────────────────────────────────
mrs_gse51092 <- msviper(signature, regulon_gse51092, nullmodel)
print(summary(mrs_gse51092))
mrs_gse51092 <- ledge(mrs_gse51092)

# ── 7. Resultados MRA ─────────────────────────────────────────────────────────
top_mrs <- viper_mrsTopTable(mrs_gse51092, n = 100)
write_tsv(top_mrs, "top100_MRs_GSE51092.txt")

pdf("Top10_MRs_GSE51092.pdf", width = 6, height = 7)
plot(mrs_gse51092, mrs = 10, cex = 0.7)
dev.off()

# ── 8. Expresión diferencial (limma) ──────────────────────────────────────────
disease_f      <- factor(disease, levels = c("control", "pSS"))
design         <- model.matrix(~ disease_f)
fit            <- lmFit(eset, design)
fit            <- eBayes(fit)
de_gse51092   <- topTable(fit, coef = 2, number = Inf,
                           adjust.method = "BH", sort.by = "P")
de_gse51092   <- tibble::rownames_to_column(de_gse51092, var = "gene")
write_tsv(de_gse51092, "DE_limma_GSE51092.txt")
cat("DE: ", sum(de_gse51092$adj.P.Val < 0.05), "genes significativos (FDR < 0.05)\n")

save(mrs_gse51092, top_mrs, regulon_gse51092, eset, disease, de_gse51092,
     file = "mra_GSE51092.rda")

packinfo(save = TRUE)
cat("Done. Resultados en top100_MRs_GSE51092.txt y DE_limma_GSE51092.txt\n")
