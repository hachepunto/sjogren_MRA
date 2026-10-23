library(RobustRankAggreg)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# Carpetas de salida
invisible(lapply(c("results/integration"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# Integración de expresión diferencial restringida a los genes medidos en AMBAS
# matrices (8,007). Sustituye la RRA sobre las tablas completas, en la que un
# gen medido en una sola plataforma podía entrar como "robusto".
#
# Tres criterios por gen (todos exigen el mismo signo de logFC en ambas cohortes):
#   (a) Replicación: FDR < 0.05 en cada cohorte por separado (criterio primario).
#   (b) Meta-análisis Stouffer sobre z derivados de p y signo de limma; BH sobre
#       los 8,007 genes.
#   (c) RRA restringida, N = genes concordantes en esa dirección; BH sobre Score.

# ── 1. Tabla de genes compartidos (de compare_ED.R) ──────────────────────────
d <- read.delim("results/de/DE_comparacion.txt")
d$same_sign <- sign(d$logFC_51) == sign(d$logFC_84)
d$direction <- ifelse(!d$same_sign, "discordant",
                      ifelse(d$logFC_51 > 0, "up", "down"))

# ── 2. Replicación por FDR en ambas cohortes ─────────────────────────────────
d$replicated <- d$same_sign & d$adj.P.Val_51 < 0.05 & d$adj.P.Val_84 < 0.05

# ── 3. Stouffer ───────────────────────────────────────────────────────────────
z_from <- function(p, lfc) sign(lfc) * qnorm(p / 2, lower.tail = FALSE)
d$z_51         <- z_from(d$P.Value_51, d$logFC_51)
d$z_84         <- z_from(d$P.Value_84, d$logFC_84)
d$z_stouffer   <- (d$z_51 + d$z_84) / sqrt(2)
d$p_stouffer   <- 2 * pnorm(-abs(d$z_stouffer))
d$BH_stouffer  <- p.adjust(d$p_stouffer, method = "BH")
d$meta_sig     <- d$same_sign & d$BH_stouffer < 0.05

# ── 4. RRA restringida a genes concordantes, por dirección ───────────────────
run_rra <- function(dd) {
    r <- aggregateRanks(list(GSE51092 = dd$gene[order(dd$P.Value_51)],
                             GSE84844 = dd$gene[order(dd$P.Value_84)]),
                        N = nrow(dd))
    r$RRA_BH <- p.adjust(r$Score, method = "BH")
    r
}
rra_up <- run_rra(d[d$direction == "up", ])
rra_dn <- run_rra(d[d$direction == "down", ])
rra_all <- rbind(rra_up, rra_dn)
d$RRA_score <- rra_all$Score[match(d$gene, rra_all$Name)]
d$RRA_BH    <- rra_all$RRA_BH[match(d$gene, rra_all$Name)]

# ── 5. Referencia: listas RRA de la v2 corregida (sobre tablas completas) ────
v2_up <- read.delim("results/de/DE_RRA_upregulated.txt")$Name
v2_dn <- read.delim("results/de/DE_RRA_downregulated.txt")$Name
d$v2_rra <- ifelse(d$gene %in% v2_up, "up", ifelse(d$gene %in% v2_dn, "down", ""))

# ── 6. Salida ─────────────────────────────────────────────────────────────────
d <- d[order(d$p_stouffer), ]
out_cols <- c("gene", "logFC_51", "P.Value_51", "adj.P.Val_51",
              "logFC_84", "P.Value_84", "adj.P.Val_84",
              "direction", "replicated", "z_stouffer", "p_stouffer",
              "BH_stouffer", "meta_sig", "RRA_score", "RRA_BH", "v2_rra")
write.table(d[, out_cols], "results/integration/DE_meta_shared.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(d, rra_up, rra_dn, file = "results/integration/meta_DE_results.rda")
