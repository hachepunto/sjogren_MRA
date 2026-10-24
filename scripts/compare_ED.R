library(RobustRankAggreg)
library(readr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# Carpetas de salida
invisible(lapply(c("results/de", "results/enrichment"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ── 1. Cargar resultados DE ───────────────────────────────────────────────────
de51 <- as.data.frame(read_tsv("results/de/DE_limma_GSE51092.txt"))
de84 <- as.data.frame(read_tsv("results/de/DE_limma_GSE84844.txt"))

# ── 2. Genes compartidos y tabla comparativa ──────────────────────────────────
shared <- merge(
    de51[, c("gene", "logFC", "P.Value", "adj.P.Val")],
    de84[, c("gene", "logFC", "P.Value", "adj.P.Val")],
    by = "gene", suffixes = c("_51", "_84")
)
# Clasificar concordancia (significativo = FDR < 0.05 en ambos, misma dirección)
shared$sig51 <- shared$adj.P.Val_51 < 0.05
shared$sig84 <- shared$adj.P.Val_84 < 0.05
shared$concordant_up   <- shared$sig51 & shared$sig84 & shared$logFC_51 > 0 & shared$logFC_84 > 0
shared$concordant_down <- shared$sig51 & shared$sig84 & shared$logFC_51 < 0 & shared$logFC_84 < 0

write.table(shared, file = "results/de/DE_comparacion.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# ── 3. Robust Rank Aggregation ────────────────────────────────────────────────
# Genes ordenados por P.Value ascendente dentro de cada dirección: en ambas
# direcciones el gen más significativo ocupa el rango 1, que es el "mejor"
# para aggregateRanks. (Versión anterior invertía la lista down con rev(),
# lo que premiaba a los genes MENOS significativos; corregido en la revisión.)

run_rra_de <- function(d1, d2, direction = c("up", "down")) {
    direction <- match.arg(direction)
    if (direction == "up") {
        l1 <- d1$gene[d1$logFC > 0][order(d1$P.Value[d1$logFC > 0])]
        l2 <- d2$gene[d2$logFC > 0][order(d2$P.Value[d2$logFC > 0])]
    } else {
        l1 <- d1$gene[d1$logFC < 0][order(d1$P.Value[d1$logFC < 0])]
        l2 <- d2$gene[d2$logFC < 0][order(d2$P.Value[d2$logFC < 0])]
    }
    rra <- aggregateRanks(list(GSE51092 = l1, GSE84844 = l2))
    rra[order(rra$Score), ]
}

rra_up <- run_rra_de(de51, de84, "up")
rra_dn <- run_rra_de(de51, de84, "down")

write.table(rra_up[rra_up$Score < 0.05, ], file = "results/de/DE_RRA_upregulated.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rra_dn[rra_dn$Score < 0.05, ], file = "results/de/DE_RRA_downregulated.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# ── 4. Genes replicados (usados en ORA y figura) ─────────────────────────────
# Criterio primario: FDR < 0.05 en ambas cohortes con el mismo signo
# (concordant_up / concordant_down de la sección 2). La RRA de la sección 3
# se conserva solo como referencia (ver meta_DE.R para Stouffer y RRA
# restringida sobre los genes compartidos).
sig_up <- shared$gene[shared$concordant_up]
sig_dn <- shared$gene[shared$concordant_down]

# z combinado (Stouffer) para ordenar los paneles C y D
z_from <- function(p, lfc) sign(lfc) * qnorm(p / 2, lower.tail = FALSE)
shared$z_stouffer <- (z_from(shared$P.Value_51, shared$logFC_51) +
                      z_from(shared$P.Value_84, shared$logFC_84)) / sqrt(2)
shared$p_stouffer <- 2 * pnorm(-abs(shared$z_stouffer))

# ── 5. ORA sobre genes DE robustos (set único, sin dirección) ────────────────
sym2entrez <- function(genes) {
    bitr(as.character(genes), fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
}

sig_all         <- union(as.character(sig_up), as.character(sig_dn))
sig_all_entrez  <- sym2entrez(sig_all)
universe_entrez <- sym2entrez(shared$gene)

ego_de <- enrichGO(
    gene          = sig_all_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
)
ego_de <- simplify(ego_de, cutoff = 0.6, by = "p.adjust")

ereact_de <- enrichPathway(
    gene          = sig_all_entrez,
    universe      = universe_entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
)

write.table(as.data.frame(ego_de),    "results/enrichment/DE_ORA_GO_BP.txt",    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(ereact_de), "results/enrichment/DE_ORA_Reactome.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ── 6. Guardar resultados ─────────────────────────────────────────────────────
write.table(shared[shared$concordant_up,   c("gene", "logFC_51", "adj.P.Val_51", "logFC_84", "adj.P.Val_84")],
            "results/de/DE_replicated_up.txt",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(shared[shared$concordant_down, c("gene", "logFC_51", "adj.P.Val_51", "logFC_84", "adj.P.Val_84")],
            "results/de/DE_replicated_down.txt", sep = "\t", quote = FALSE, row.names = FALSE)
save(shared, rra_up, rra_dn, sig_up, sig_dn, ego_de, ereact_de,
     file = "results/de/comparison_DE_results.rda")
