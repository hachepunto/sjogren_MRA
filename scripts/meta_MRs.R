library(RobustRankAggreg)

# Carpetas de salida
invisible(lapply(c("results/integration"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# Integración de MRs restringida al universo compartido (TFs con regulón en
# AMBAS redes). Sustituye la RRA sobre la unión de reguladores, en la que un TF
# presente en una sola red podía alcanzar Score < 0.05 por sí solo.
#
# Dos métricas por TF:
#   (a) Meta-análisis de actividad (Stouffer): los NES de msviper son z-scores
#       bajo el modelo nulo, así que z = (NES1 + NES2)/sqrt(2); BH sobre los
#       240 TFs compartidos. Se exige además el mismo signo en ambas cohortes.
#   (b) RRA restringida: listas solo con TFs de signo concordante, N = número
#       de TFs en esa dirección, BH sobre el Score.

# ── 1. Tabla de TFs compartidos (de compare_MRs.R) ───────────────────────────
sh <- read.delim("results/integration/MRs_comparacion.txt")
n1 <- "NES_gse51092"; n2 <- "NES_84844"

sh$same_sign <- sign(sh[[n1]]) == sign(sh[[n2]])
sh$direction <- ifelse(!sh$same_sign, "discordant",
                       ifelse(sh[[n1]] > 0, "activated", "repressed"))

# ── 2. Stouffer sobre NES ─────────────────────────────────────────────────────
sh$z_stouffer  <- (sh[[n1]] + sh[[n2]]) / sqrt(2)
sh$p_stouffer  <- 2 * pnorm(-abs(sh$z_stouffer))
sh$BH_stouffer <- p.adjust(sh$p_stouffer, method = "BH")

# ── 3. RRA restringida a TFs concordantes, por dirección ─────────────────────
run_rra <- function(d, decreasing) {
    o1 <- order(d[[n1]], decreasing = decreasing)
    o2 <- order(d[[n2]], decreasing = decreasing)
    r  <- aggregateRanks(list(GSE51092 = d$Regulon[o1], GSE84844 = d$Regulon[o2]),
                         N = nrow(d))
    r$RRA_BH <- p.adjust(r$Score, method = "BH")
    r
}
rra_act <- run_rra(sh[sh$direction == "activated", ], decreasing = TRUE)
rra_rep <- run_rra(sh[sh$direction == "repressed", ], decreasing = FALSE)
rra_all <- rbind(rra_act, rra_rep)
sh$RRA_score <- rra_all$Score[match(sh$Regulon, rra_all$Name)]
sh$RRA_BH    <- rra_all$RRA_BH[match(sh$Regulon, rra_all$Name)]

# ── 4. Referencia: lista de consenso de la v2 (RRA sobre la unión) ───────────
v2_act <- read.delim("results/integration/MRs_RRA_activated.txt")$Name
v2_rep <- read.delim("results/integration/MRs_RRA_repressed.txt")$Name
sh$v2_consensus <- ifelse(sh$Regulon %in% v2_act, "activated",
                   ifelse(sh$Regulon %in% v2_rep, "repressed", ""))

# ── 5. Concordancia global entre cohortes ────────────────────────────────────
rho <- cor.test(sh[[n1]], sh[[n2]], method = "spearman", exact = FALSE)
# Enriquecimiento del top-25 de GSE51092 en el ranking de GSE84844 (y viceversa):
# Wilcoxon rank-sum de los NES de la otra cohorte, top-25 vs resto.
top_test <- function(ref, other, k = 25) {
    top <- order(sh[[ref]], decreasing = TRUE)[seq_len(k)]
    wilcox.test(sh[[other]][top], sh[[other]][-top], alternative = "greater")
}
w_51_in_84 <- top_test(n1, n2)
w_84_in_51 <- top_test(n2, n1)
concordance <- data.frame(
    test      = c("Spearman rho (240 TFs)",
                  "Top-25 GSE51092 in GSE84844 (Wilcoxon, greater)",
                  "Top-25 GSE84844 in GSE51092 (Wilcoxon, greater)"),
    statistic = c(rho$estimate, w_51_in_84$statistic, w_84_in_51$statistic),
    p.value   = c(rho$p.value, w_51_in_84$p.value, w_84_in_51$p.value)
)

# ── 6. Salida ─────────────────────────────────────────────────────────────────
sh <- sh[order(sh$p_stouffer), ]
out_cols <- c("Regulon", n1, "p.value_gse51092", "rank_gse51092",
              n2, "p.value_84844", "rank_84844",
              "direction", "z_stouffer", "p_stouffer", "BH_stouffer",
              "RRA_score", "RRA_BH", "v2_consensus")
write.table(sh[, out_cols], "results/integration/MRs_meta_shared.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concordance, "results/integration/MRs_concordance_tests.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(sh, concordance, rra_act, rra_rep, file = "results/integration/meta_MRs_results.rda")
