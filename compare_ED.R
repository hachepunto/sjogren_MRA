library(RobustRankAggreg)
library(readr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# ── 1. Cargar resultados DE ───────────────────────────────────────────────────
de51 <- as.data.frame(read_tsv("DE_limma_GSE51092.txt"))
de84 <- as.data.frame(read_tsv("DE_limma_GSE84844.txt"))


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

write.table(shared, file = "DE_comparacion.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# ── 3. Robust Rank Aggregation ────────────────────────────────────────────────
# Genes ordenados por P.Value ascendente dentro de cada dirección.
# Para down se invierte la lista (mismo criterio que compare_MRs.R):
# los genes con logFC < 0 aparecen en orden P.Value ascendente al filtrar,
# rev() pone el más significativo al final → rank N (el "mejor" para RRA repressed).

run_rra_de <- function(d1, d2, direction = c("up", "down")) {
    direction <- match.arg(direction)
    if (direction == "up") {
        l1 <- d1$gene[d1$logFC > 0][order(d1$P.Value[d1$logFC > 0])]
        l2 <- d2$gene[d2$logFC > 0][order(d2$P.Value[d2$logFC > 0])]
    } else {
        l1 <- rev(d1$gene[d1$logFC < 0][order(d1$P.Value[d1$logFC < 0])])
        l2 <- rev(d2$gene[d2$logFC < 0][order(d2$P.Value[d2$logFC < 0])])
    }
    rra <- aggregateRanks(list(GSE51092 = l1, GSE84844 = l2))
    rra[order(rra$Score), ]
}

rra_up <- run_rra_de(de51, de84, "up")
rra_dn <- run_rra_de(de51, de84, "down")

write.table(rra_up[rra_up$Score < 0.05, ], file = "DE_RRA_upregulated.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rra_dn[rra_dn$Score < 0.05, ], file = "DE_RRA_downregulated.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# ── 4. Genes robustos DE (usados en ORA y figura) ────────────────────────────
sig_up <- rra_up$Name[rra_up$Score < 0.05]
sig_dn <- rra_dn$Name[rra_dn$Score < 0.05]

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

write.table(as.data.frame(ego_de),    "DE_ORA_GO_BP.txt",    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(ereact_de), "DE_ORA_Reactome.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ── 6. Guardar resultados ─────────────────────────────────────────────────────
save(shared, rra_up, rra_dn, ego_de, ereact_de,
     file = "comparison_DE_results.rda")

# ── 7. Figure (ggplot2 + patchwork) ──────────────────────────────────────────
# Layout (column-major, nrow=2):   A | C
#                                  B | D
# A = scatter (1:1), B = ORA (1:1), C = Up lollipop (0.5:1), D = Down lollipop (0.5:1)
# Column widths: left=1, right=0.5 → PDF 12×16 gives 8×8 left panels, 4×8 right panels

shared$status <- "Not significant"
shared$status[(shared$sig51 | shared$sig84) & !shared$concordant_up & !shared$concordant_down] <- "Single dataset only"
shared$status[shared$gene %in% sig_up] <- "Up robust (RRA)"
shared$status[shared$gene %in% sig_dn] <- "Down robust (RRA)"
shared$status <- factor(shared$status,
    levels = c("Up robust (RRA)", "Down robust (RRA)", "Single dataset only", "Not significant"))

cols_status <- c(
    "Up robust (RRA)"     = "#D73027",
    "Down robust (RRA)"   = "#4575B4",
    "Single dataset only" = "#FEE090",
    "Not significant"     = "grey80"
)
shared_sig <- shared[shared$gene %in% c(sig_up, sig_dn), ]

# Panel A — logFC scatter (1:1, coord_fixed forces equal axis scales)
pA <- ggplot(shared, aes(logFC_51, logFC_84, fill = status)) +
    geom_point(shape = 21, color = "white", size = 1.2, alpha = 0.7) +
    geom_text(data = shared_sig, aes(label = gene, color = status),
              size = 1.8, vjust = -0.4, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    scale_fill_manual(values = cols_status, name = NULL) +
    scale_color_manual(values = cols_status) +
    coord_fixed() +
    labs(x = "logFC — GSE51092 (n=222)", y = "logFC — GSE84844 (n=60)",
         title = "A   Cross-dataset differential expression") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          plot.title      = element_text(face = "bold", size = 10))

# Panel B — ORA combined lollipop with bubbles (1:1)
ora_go             <- as.data.frame(ego_de)[,    c("Description", "GeneRatio", "Count", "p.adjust")]
ora_go$Database    <- "GO BP"
ora_react          <- as.data.frame(ereact_de)[, c("Description", "GeneRatio", "Count", "p.adjust")]
ora_react$Database <- "Reactome"
ora_df             <- rbind(ora_go, ora_react)
ora_df$GeneRatioNum <- sapply(strsplit(ora_df$GeneRatio, "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]))
ora_df$Description  <- factor(ora_df$Description,
                               levels = ora_df$Description[order(ora_df$GeneRatioNum)])

pB <- ggplot(ora_df, aes(x = GeneRatioNum, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatioNum, yend = Description, color = Database),
                 linewidth = 1.2) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = "white") +
    scale_color_manual(values = c("GO BP" = "#E08533", "Reactome" = "#7B3F9E"),
                       name = "Database") +
    scale_fill_gradient(low = "#4575B4", high = "#FEE090", name = "FDR") +
    scale_size_continuous(name = "Count", range = c(3, 8)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    guides(
        color = guide_legend(order = 1),
        fill  = guide_colorbar(order = 2),
        size  = guide_legend(order = 3,
                             override.aes = list(fill = "grey50", color = "white"))
    ) +
    labs(x = "Gene Ratio", y = NULL,
         title = "B   Functional enrichment — robust DE genes") +
    theme_bw(base_size = 9) +
    theme(axis.text.y    = element_text(size = 8),
          plot.title     = element_text(face = "bold", size = 10),
          legend.position = "left",
          aspect.ratio    = 1)

# Panel C — Up-regulated lollipop (0.5:1)
lol_up      <- head(rra_up[rra_up$Score < 0.05, ], 30)
lol_up$logp <- -log10(lol_up$Score)

pC <- ggplot(lol_up, aes(x = logp, y = reorder(Name, logp))) +
    geom_segment(aes(x = 0, xend = logp, yend = reorder(Name, logp)),
                 color = "#D73027", linewidth = 1) +
    geom_point(fill = "#D73027", color = "white", shape = 21, size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = expression(-log[10](Score)), y = NULL,
         title = "C   Up-regulated in pSS (RRA)") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7),
          plot.title  = element_text(face = "bold", size = 10))

# Panel D — Down-regulated lollipop (0.5:1)
lol_dn      <- head(rra_dn[rra_dn$Score < 0.05, ], 30)
lol_dn$logp <- -log10(lol_dn$Score)

pD <- ggplot(lol_dn, aes(x = logp, y = reorder(Name, logp))) +
    geom_segment(aes(x = 0, xend = logp, yend = reorder(Name, logp)),
                 color = "#4575B4", linewidth = 1) +
    geom_point(fill = "#4575B4", color = "white", shape = 21, size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = expression(-log[10](Score)), y = NULL,
         title = "D   Down-regulated in pSS (RRA)") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7),
          plot.title  = element_text(face = "bold", size = 10))

fig <- pA + pC + pB + pD +
    plot_layout(ncol = 2, widths = c(1, 0.5), heights = c(1, 1))

pdf("DE_integrative_figure.pdf", width = 16, height = 12)
print(fig)
dev.off()


