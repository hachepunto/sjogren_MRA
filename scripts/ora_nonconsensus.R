library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# ORA de los regulones ARACNE-AP de los cinco reguladores que no entran al
# consenso pero se evalúan aparte con la red fija (Métodos Suplementarios 2.1).
# Mismo fondo, umbrales y reducción semántica que ora_regulons.R.
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

load("results/mra/mra_GSE51092.rda")                 # regulon_gse51092, eset

MRS <- c("HESX1", "ZNF684", "THRA", "MEF2D", "FOXJ2")
stopifnot(all(MRS %in% names(regulon_gse51092)))

gene_lists_sym <- lapply(setNames(MRS, MRS),
                         function(tf) names(regulon_gse51092[[tf]]$tfmode))
cat("Tamaño de regulones:\n")
print(sapply(gene_lists_sym, length))

# Genes de reparación y desubiquitinación citados en la Sección 4.3 para FOXJ2
foxj2_cited <- c("SMARCA4", "CHD4", "RPA1", "XRCC6", "REV3L", "VCP", "USP7")
cat("Genes citados en el regulón de FOXJ2:\n")
print(setNames(foxj2_cited %in% gene_lists_sym$FOXJ2, foxj2_cited))

sym2entrez <- function(genes) {
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
}
gene_lists_entrez <- lapply(gene_lists_sym, sym2entrez)
universe_entrez   <- sym2entrez(rownames(eset))

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

# Una sola tabla larga; compareCluster devuelve NULL si ningún regulón enriquece
as_tab <- function(cc, db) {
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) return(NULL)
    x <- as.data.frame(cc)
    data.frame(database = db, x[, c("Cluster", "ID", "Description", "GeneRatio", "BgRatio",
                                    "pvalue", "p.adjust", "qvalue", "geneID", "Count")])
}
ora <- rbind(as_tab(cc_go, "GO BP"), as_tab(cc_react, "Reactome"))
if (is.null(ora)) ora <- data.frame(database = character(), Cluster = character())
names(ora)[names(ora) == "Cluster"] <- "regulator"
cat("Términos por regulador y base:\n")
print(table(factor(ora$regulator, levels = MRS), ora$database))

write.table(ora, "results/enrichment/ORA_nonconsensus.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
save(cc_go, cc_react, gene_lists_sym, file = "results/enrichment/ORA_nonconsensus.rda")
