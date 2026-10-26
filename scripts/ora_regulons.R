library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

source("scripts/helpers.R")

# Carpetas de salida
invisible(lapply(c("results/enrichment", "results/figures"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

load("results/mra/mra_GSE51092.rda")                 # regulon_gse51092, eset

# ── 1. MRs a analizar: consenso entre cohortes (meta_MRs.R) ──────────────────
# Regulones ARACNe-AP de GSE51092; fondo = genes de la matriz de expresión.
cons    <- read_consensus()
mrs_act <- cons$Regulon[cons$direction == "activated"]
mrs_rep <- cons$Regulon[cons$direction == "repressed"]
all_mrs <- c(mrs_act, mrs_rep)

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

# "\u00ad" (guion suave): pdf() dibuja "-" como signo menos
con_guion <- function(cc) {
    cc@compareClusterResult$Description <- gsub("-", "\u00ad", cc@compareClusterResult$Description, fixed = TRUE)
    cc
}

p_go <- dotplot(con_guion(cc_go), showCategory = 5, font.size = 9) +
    labs(title = "A   GO Biological Process") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold", size = 10))

p_react <- dotplot(con_guion(cc_react), showCategory = 5, font.size = 9) +
    labs(title = "B   Reactome") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold", size = 10))

combined <- p_go | p_react + plot_layout(widths = c(3, 1))

pdf("results/figures/ORA_combined.pdf", width = 16, height = 8)
print(combined)
dev.off()

# ── 7. Guardar resultados ─────────────────────────────────────────────────────
save(cc_go, cc_react, file = "results/enrichment/ORA_results.rda")
write.table(as.data.frame(cc_go),    "results/enrichment/ORA_GO_BP.txt",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(cc_react), "results/enrichment/ORA_Reactome.txt", sep = "\t", quote = FALSE, row.names = FALSE)

