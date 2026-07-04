library(GEOquery)
library(matrixStats)
library(limma)
library(tibble)
library(readr)
library(openxlsx)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# ── 1. Expression matrix (47323 probes → genes × 160 samples) ────────────────
gse     <- getGEO(filename = "data/GSE66795_series_matrix.txt", GSEMatrix = TRUE)
eset    <- exprs(gse)
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

# ── 2. Disease labels ─────────────────────────────────────────────────────────
# pData[["patient group:ch1"]]: "Control" (29) | "Patient" (131, all fatigue subgroups)
patient_group <- pData(gse)[["patient group:ch1"]]
disease       <- ifelse(patient_group == "Control", "control", "pSS")
print(table(disease))

# ── 3. Differential expression (limma/eBayes) ─────────────────────────────────
disease_f <- factor(disease, levels = c("control", "pSS"))
design    <- model.matrix(~ disease_f)
fit       <- lmFit(eset, design)
fit       <- eBayes(fit)
de_gse66795 <- topTable(fit, coef = 2, number = Inf,
                         adjust.method = "BH", sort.by = "P")
de_gse66795 <- rownames_to_column(de_gse66795, var = "gene")
write_tsv(de_gse66795, "DE_limma_GSE66795.txt")

# ── 4. Overlap with RRA-robust DE gene sets ───────────────────────────────────
rra_up <- read_tsv("DE_RRA_upregulated.txt",   show_col_types = FALSE)
rra_dn <- read_tsv("DE_RRA_downregulated.txt",  show_col_types = FALSE)

sig_up66 <- de_gse66795$gene[de_gse66795$adj.P.Val < 0.05 & de_gse66795$logFC > 0]
sig_dn66 <- de_gse66795$gene[de_gse66795$adj.P.Val < 0.05 & de_gse66795$logFC < 0]

overlap_up <- rra_up[rra_up$Name %in% sig_up66, ]
overlap_dn <- rra_dn[rra_dn$Name %in% sig_dn66, ]

print(data.frame(
    set        = c("RRA_up (n=636)", "RRA_down (n=514)"),
    replicated = c(nrow(overlap_up), nrow(overlap_dn)),
    pct        = round(c(nrow(overlap_up) / nrow(rra_up),
                         nrow(overlap_dn) / nrow(rra_dn)) * 100, 1)
))

# ── 5. Save ───────────────────────────────────────────────────────────────────
build_table_s1 <- function(overlap, de) {
    out <- merge(overlap, de[, c("gene", "logFC", "adj.P.Val")],
                 by.x = "Name", by.y = "gene", all.x = TRUE, sort = FALSE)
    names(out)[names(out) == "Score"] <- "RRA_score"
    out <- out[order(out$RRA_score), c("Name", "RRA_score", "logFC", "adj.P.Val")]
    out
}

table_s1_up <- build_table_s1(overlap_up, de_gse66795)
table_s1_dn <- build_table_s1(overlap_dn, de_gse66795)

wb <- createWorkbook()

addWorksheet(wb, "Upregulated_165")
writeData(wb, "Upregulated_165",
          x = "Supplementary Table S1a. Genes from the cross-cohort IFN-upregulated RRA signature (n=636) that replicate in GSE66795 (limma/eBayes, adj.P.Val < 0.05 and logFC > 0; Illumina HumanHT-12 v4.0).",
          startRow = 1, startCol = 1)
writeData(wb, "Upregulated_165", table_s1_up, startRow = 2)

addWorksheet(wb, "Downregulated_11")
writeData(wb, "Downregulated_11",
          x = "Supplementary Table S1b. Genes from the cross-cohort IFN-downregulated RRA signature (n=514) that replicate in GSE66795 (limma/eBayes, adj.P.Val < 0.05 and logFC < 0; Illumina HumanHT-12 v4.0).",
          startRow = 1, startCol = 1)
writeData(wb, "Downregulated_11", table_s1_dn, startRow = 2)

saveWorkbook(wb, "Table_S1_DE_GSE66795_replication.xlsx", overwrite = TRUE)

save(eset, disease, de_gse66795, overlap_up, overlap_dn,
     file = "DE_GSE66795_replication.rda")
