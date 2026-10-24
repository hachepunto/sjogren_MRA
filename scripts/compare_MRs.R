library(viper)
library(RobustRankAggreg)

# Carpetas de salida
invisible(lapply(c("results/integration", "results/figures"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

load("results/mra/mra_GSE51092.rda")   # mrs_gse51092, regulon_gse51092, eset, disease
load("results/mra/mra_GSE84844.rda")   # mrs_gse84844

# ── 1. Extraer resultados completos de ambos objetos msviper ──────────────────
extract_mrs <- function(mrs) {
    nes  <- mrs$es$nes
    pval <- mrs$es$p.value
    data.frame(
        Regulon = names(nes),
        NES     = nes,
        p.value = pval,
        row.names = NULL
    ) |> (\(d) d[order(-d$NES), ])()
}

sum_gse51092 <- extract_mrs(mrs_gse51092)
sum_84844    <- extract_mrs(mrs_gse84844)

sum_gse51092$rank_gse51092 <- seq_len(nrow(sum_gse51092))
sum_84844$rank_84844       <- seq_len(nrow(sum_84844))

# ── 2. Tabla comparativa de MRs compartidos ───────────────────────────────────
shared <- merge(
    sum_gse51092[, c("Regulon", "NES", "p.value", "rank_gse51092")],
    sum_84844[,   c("Regulon", "NES", "p.value", "rank_84844")],
    by = "Regulon",
    suffixes = c("_gse51092", "_84844")
)
shared <- shared[order(shared$rank_gse51092), ]

cat("MRs en GSE51092:            ", nrow(sum_gse51092), "\n")
cat("MRs en GSE84844:            ", nrow(sum_84844), "\n")
cat("MRs compartidos:            ", nrow(shared), "\n")
cat("Compartidos top-25 en ambos:",
    sum(shared$rank_gse51092 <= 25 & shared$rank_84844 <= 25), "\n\n")

# Top compartidos
print(head(shared[shared$rank_gse51092 <= 25 & shared$rank_84844 <= 25, ], 20))

write.table(shared, file = "results/integration/MRs_comparacion.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# ── 3. Scatter NES: GSE51092 vs GSE84844 ─────────────────────────────────────
top_both <- shared$rank_gse51092 <= 25 | shared$rank_84844 <= 25

pdf("results/figures/MRs_scatter_NES.pdf", width = 6, height = 6)
plot(shared$NES_gse51092, shared$NES_84844,
     pch  = 16,
     col  = ifelse(shared$rank_gse51092 <= 25 & shared$rank_84844 <= 25,
                   "firebrick", "grey70"),
     xlab = "NES — GSE51092 (ARACNE-AP)",
     ylab = "NES — GSE84844 (corto)",
     main = "Master Regulators: GSE51092 vs GSE84844")
abline(h = 0, v = 0, lty = 2, col = "grey50")
# Etiquetar MRs top en ambos datasets
top_label <- shared[shared$rank_gse51092 <= 25 & shared$rank_84844 <= 25, ]
text(top_label$NES_gse51092, top_label$NES_84844,
     labels = top_label$Regulon, cex = 0.7, pos = 3)
dev.off()

# ── 4. Diagrama de Venn (top 25 de cada dataset) ─────────────────────────────
top25_gse51092 <- sum_gse51092$Regulon[1:25]
top25_84844    <- sum_84844$Regulon[1:25]
only_gse51092  <- setdiff(top25_gse51092, top25_84844)
only_84844     <- setdiff(top25_84844, top25_gse51092)
in_both        <- intersect(top25_gse51092, top25_84844)

cat("\nTop-25 solo en GSE51092: ", paste(only_gse51092, collapse = ", "), "\n")
cat("Top-25 solo en GSE84844: ", paste(only_84844,     collapse = ", "), "\n")
cat("Top-25 en ambos:         ", paste(in_both,         collapse = ", "), "\n")

# ── 5. Robust Rank Aggregation ────────────────────────────────────────────────
run_rra <- function(s1, s2, direction = c("activated", "repressed")) {
    direction <- match.arg(direction)
    if (direction == "activated") {
        l1 <- s1$Regulon[s1$NES > 0]
        l2 <- s2$Regulon[s2$NES > 0]
    } else {
        l1 <- rev(s1$Regulon[s1$NES < 0])
        l2 <- rev(s2$Regulon[s2$NES < 0])
    }
    rra <- aggregateRanks(list(GSE51092 = l1, GSE84844 = l2))
    rra[order(rra$Score), ]
}

rra_act <- run_rra(sum_gse51092, sum_84844, "activated")
rra_rep <- run_rra(sum_gse51092, sum_84844, "repressed")

cat("\n── RRA MRs activados en pSS (Score < 0.05) ──\n")
print(head(rra_act[rra_act$Score < 0.05, ], 30))

cat("\n── RRA MRs reprimidos en pSS (Score < 0.05) ──\n")
print(head(rra_rep[rra_rep$Score < 0.05, ], 30))

write.table(rra_act[rra_act$Score < 0.05, ], file = "results/integration/MRs_RRA_activated.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rra_rep[rra_rep$Score < 0.05, ], file = "results/integration/MRs_RRA_repressed.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# ── 6. Guardar resultados ─────────────────────────────────────────────────────
save(sum_gse51092, sum_84844, shared, rra_act, rra_rep,
     file = "results/integration/comparison_results.rda")
