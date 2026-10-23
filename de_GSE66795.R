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

# ── 4. Solapamiento con los conjuntos replicados (FDR < 0.05 en ambas cohortes)
# Se conservan los nombres rra_up / rra_dn por compatibilidad con el resto del
# script; la columna Name es el símbolo y Score el p combinado (Stouffer).
load("comparison_DE_results.rda")   # shared (con z_stouffer, p_stouffer), sig_up, sig_dn
rep_cols <- c("gene", "p_stouffer", "logFC_51", "adj.P.Val_51", "logFC_84", "adj.P.Val_84")
rra_up <- shared[shared$concordant_up,   rep_cols]
rra_dn <- shared[shared$concordant_down, rep_cols]
names(rra_up)[1:2] <- names(rra_dn)[1:2] <- c("Name", "Score")

sig_up66 <- de_gse66795$gene[de_gse66795$adj.P.Val < 0.05 & de_gse66795$logFC > 0]
sig_dn66 <- de_gse66795$gene[de_gse66795$adj.P.Val < 0.05 & de_gse66795$logFC < 0]

overlap_up <- rra_up[rra_up$Name %in% sig_up66, ]
overlap_dn <- rra_dn[rra_dn$Name %in% sig_dn66, ]

print(data.frame(
    set        = c(sprintf("replicated_up (n=%d)", nrow(rra_up)), sprintf("replicated_down (n=%d)", nrow(rra_dn))),
    replicated = c(nrow(overlap_up), nrow(overlap_dn)),
    pct        = round(c(nrow(overlap_up) / nrow(rra_up),
                         nrow(overlap_dn) / nrow(rra_dn)) * 100, 1)
))

# ── 5. Save ───────────────────────────────────────────────────────────────────
build_table_s1 <- function(overlap, de) {
    out <- merge(overlap, de[, c("gene", "logFC", "adj.P.Val")],
                 by.x = "Name", by.y = "gene", all.x = TRUE, sort = FALSE)
    names(out)[names(out) == "Score"] <- "p_stouffer_discovery_validation"
    names(out)[names(out) == "logFC"] <- "logFC_GSE66795"
    names(out)[names(out) == "adj.P.Val"] <- "adj.P.Val_GSE66795"
    out <- out[order(out$p_stouffer_discovery_validation),
               c("Name", "logFC_51", "adj.P.Val_51", "logFC_84", "adj.P.Val_84",
                 "p_stouffer_discovery_validation", "logFC_GSE66795", "adj.P.Val_GSE66795")]
    out
}

table_s1_up <- build_table_s1(overlap_up, de_gse66795)
table_s1_dn <- build_table_s1(overlap_dn, de_gse66795)

wb <- createWorkbook()

sheet_up <- sprintf("Upregulated_%d",   nrow(table_s1_up))
sheet_dn <- sprintf("Downregulated_%d", nrow(table_s1_dn))

addWorksheet(wb, sheet_up)
writeData(wb, sheet_up,
          x = sprintf("Supplementary Table S1a. Upregulated genes replicated in GSE51092 and GSE84844 (FDR < 0.05 in both, n=%d) that also replicate in GSE66795 (limma/eBayes, adj.P.Val < 0.05 and logFC > 0; Illumina HumanHT-12 v4.0).", nrow(rra_up)),
          startRow = 1, startCol = 1)
writeData(wb, sheet_up, table_s1_up, startRow = 2)

addWorksheet(wb, sheet_dn)
writeData(wb, sheet_dn,
          x = sprintf("Supplementary Table S1b. Downregulated genes replicated in GSE51092 and GSE84844 (FDR < 0.05 in both, n=%d) that also replicate in GSE66795 (limma/eBayes, adj.P.Val < 0.05 and logFC < 0; Illumina HumanHT-12 v4.0).", nrow(rra_dn)),
          startRow = 1, startCol = 1)
writeData(wb, sheet_dn, table_s1_dn, startRow = 2)

saveWorkbook(wb, "Table_S1_DE_GSE66795_replication.xlsx", overwrite = TRUE)

save(eset, disease, de_gse66795, overlap_up, overlap_dn,
     file = "DE_GSE66795_replication.rda")
