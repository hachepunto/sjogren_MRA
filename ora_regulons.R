library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

load("mra_GSE51092.rda")        # regulon_gse51092, eset
load("comparison_results.rda")  # rra_act, rra_rep

# ── 1. MRs a analizar (top 5 por dirección) ───────────────────────────────────
top10_act <- as.character(rra_act$Name[1:10])
top10_rep <- as.character(rra_rep$Name[1:10])
all_mrs  <- c(top10_act, top10_rep)

cat("Activados:", paste(top10_act, collapse = ", "), "\n")
cat("Reprimidos:", paste(top10_rep, collapse = ", "), "\n")

# Verificar que están en el regulon
missing <- all_mrs[!all_mrs %in% names(regulon_gse51092)]
if (length(missing) > 0) cat("Aviso — no en regulon:", paste(missing, collapse = ", "), "\n")
all_mrs <- all_mrs[all_mrs %in% names(regulon_gse51092)]

# ── 2. Targets de cada regulon (gene symbols) ─────────────────────────────────
gene_lists_sym <- lapply(setNames(all_mrs, all_mrs),
                         function(tf) names(regulon_gse51092[[tf]]$tfmode))

cat("Tamaño de regulones:\n")
print(sapply(gene_lists_sym, length))

# ── 3. Conversión a Entrez IDs ────────────────────────────────────────────────
sym2entrez <- function(genes) {
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
}

gene_lists_entrez <- lapply(gene_lists_sym, sym2entrez)
universe_entrez   <- sym2entrez(rownames(eset))

# ── 4. ORA — GO Biological Process ───────────────────────────────────────────
cc_go <- compareCluster(
    geneClusters  = gene_lists_entrez,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    universe      = universe_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
)
cc_go <- clusterProfiler::simplify(cc_go, cutoff = 0.6, by = "p.adjust")

# ── 5. ORA — Reactome ─────────────────────────────────────────────────────────
cc_react <- compareCluster(
    geneClusters  = gene_lists_entrez,
    fun           = "enrichPathway",
    organism      = "human",
    universe      = universe_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
)

# ── 6. Figura combinada (paneles A y B) ──────────────────────────────────────
library(ggplot2)
library(patchwork)

p_go <- dotplot(cc_go, showCategory = 5, font.size = 9) +
    labs(title = "A   GO Biological Process") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold", size = 10))

p_react <- dotplot(cc_react, showCategory = 5, font.size = 9) +
    labs(title = "B   Reactome") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold", size = 10))

combined <- p_go | p_react + plot_layout(widths = c(3, 1))

pdf("ORA_combined.pdf", width = 16, height = 8)
print(combined)
dev.off()

# ── 7. Guardar resultados ─────────────────────────────────────────────────────
save(cc_go, cc_react, file = "ORA_results.rda")
write.table(as.data.frame(cc_go),    "ORA_GO_BP.txt",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(cc_react), "ORA_Reactome.txt", sep = "\t", quote = FALSE, row.names = FALSE)

