# Libro de Tablas Suplementarias S1–S7: una hoja por tabla de results/, con su
# leyenda, más una hoja índice.
suppressMessages(library(openxlsx))
out <- "results/Supplementary_Tables_S1-S7.xlsx"

rd <- function(f) {
    x <- read.delim(f, check.names = FALSE)
    if (names(x)[1] == "") names(x)[1] <- "sample"
    x
}
# Cifras de las leyendas, calculadas a partir de las tablas
mr  <- rd("results/integration/MRs_meta_shared.txt")
de  <- rd("results/integration/DE_meta_shared.txt")
n_up <- sum(de$replicated & de$direction == "up"); n_dn <- sum(de$replicated & de$direction == "down")
fmt  <- function(x) format(x, big.mark = ",")

s1 <- "results/validation/Table_S1_DE_GSE66795_replication.xlsx"
s1x <- function(i) read.xlsx(s1, sheet = i, startRow = 2)
s1t <- function(i) sub("^Supplementary Table S1[ab]\\. ", "",   # el prefijo se antepone al escribir
                        read.xlsx(s1, sheet = i, rows = 1, colNames = FALSE)[1, 1])

sheets <- list(
  list("S1_up", "S1", s1t(1), s1x(1)),
  list("S1_down", "S1", s1t(2), s1x(2)),
  list("S2", "S2", sprintf("Cross-cohort integration of regulator activity for the %d shared transcription factors: NES, p and rank per cohort, direction, Stouffer z, p and FDR, consensus flag (same sign in both cohorts and BH_stouffer < 0.10; %d regulators), restricted RRA score and FDR (descriptive), and membership in the RRA lists computed over the union of regulators, the integration criterion of the submitted version (RRA_union_list).", nrow(mr), sum(mr$consensus)),
       mr),
  list("S3", "S3", sprintf("Cross-cohort integration of differential expression for the %s shared genes: limma statistics per cohort (_51 = GSE51092, _84 = GSE84844), direction, replication flag (FDR < 0.05 in both cohorts, same sign; %d up, %d down), Stouffer z, p and FDR, meta_sig (same sign and BH_stouffer < 0.05), restricted RRA score and FDR (descriptive), and membership in the RRA lists computed over the full tables (RRA_union_list; submitted version's criterion with the ordering of the down list corrected).", fmt(nrow(de)), n_up, n_dn),
       de),
  list("S4_regulators", "S4", "Cross-network VIPER: own-network and cross-network NES for every regulator, and fixed-network meta-analysis (aracne_*: ARACNe-AP network of GSE51092 applied to both signatures; corto_*: corto network of GSE84844 applied to both signatures).",
       rd("results/robustness/viper_cross_network.txt")),
  list("S4_concordance", "S4", "Cross-network VIPER: Spearman concordance between own-network and cross-network NES.",
       rd("results/robustness/viper_cross_network_concordance.txt")),
  list("S5_regulators", "S5", "VIPER with CollecTRI regulons (minsize 10): NES and p per cohort, regulon sizes and Stouffer meta-analysis.",
       rd("results/robustness/viper_collectri.txt")),
  list("S5_consensus", "S5", "Coverage of the 240 shared regulators in CollecTRI and their CollecTRI-based activity.",
       rd("results/robustness/consensus_vs_collectri.txt")),
  list("S5_concordance", "S5", "Spearman concordance between inferred-network and CollecTRI combined z.",
       rd("results/robustness/viper_collectri_concordance.txt")),
  list("S6_xCell_GSE51092", "S6", "xCell scores per sample, GSE51092.", rd("results/robustness/xcell_GSE51092.txt")),
  list("S6_xCell_GSE84844", "S6", "xCell scores per sample, GSE84844.", rd("results/robustness/xcell_GSE84844.txt")),
  list("S6_MCP_GSE51092", "S6", "MCP-counter scores per sample, GSE51092.", rd("results/robustness/mcpcounter_GSE51092.txt")),
  list("S6_MCP_GSE84844", "S6", "MCP-counter scores per sample, GSE84844.", rd("results/robustness/mcpcounter_GSE84844.txt")),
  list("S6_composition", "S6", "Composition tests: pSS vs control per cell type (Wilcoxon), BH across cell types within cohort.",
       rd("results/robustness/deconv_composition_tests.txt")),
  list("S6_adjusted_MRA", "S6", "Composition-adjusted versus unadjusted NES for every regulator (xCell scores of six populations regressed out of the expression, with disease status in the design so that the pSS effect is preserved), with unadjusted and adjusted Stouffer meta-analysis.",
       rd("results/robustness/mra_adjusted_composition.txt")),
  list("S6_concordance", "S6", "Spearman concordance between adjusted and unadjusted NES per cohort.",
       rd("results/robustness/mra_adjusted_concordance.txt")),
  list("S7_fgsea", "S7", "fgsea (multilevel) enrichment of positive- and negative-mode regulon targets on each cohort signature; NES_oriented flips negative-mode sets so that both reflect regulator activity.",
       rd("results/enrichment/fgsea_regulons_long.txt")),
  list("S7_concordance", "S7", "Spearman concordance between fgsea NES (positive-mode sets) and VIPER NES per cohort.",
       rd("results/enrichment/fgsea_vs_viper_concordance.txt"))
)

legends <- c(
  S1 = sprintf("Replication of the %d upregulated and %d downregulated genes in GSE66795 (%d and %d recovered; limma statistics in the three cohorts).",
               n_up, n_dn, nrow(sheets[[1]][[4]]), nrow(sheets[[2]][[4]])),
  S2 = sprintf("Cross-cohort integration of regulator activity for the %d shared transcription factors.", nrow(mr)),
  S3 = sprintf("Cross-cohort integration of differential expression for the %s shared genes.", fmt(nrow(de))),
  S4 = "Cross-network VIPER: own-network and cross-network NES for every regulator, fixed-network meta-analysis, and concordance statistics.",
  S5 = "VIPER with CollecTRI regulons: NES and p per cohort, regulon sizes, Stouffer meta-analysis, and coverage of the consensus regulators in CollecTRI.",
  S6 = "Cell-type deconvolution (xCell and MCP-counter scores per sample), composition tests, and composition-adjusted versus unadjusted NES for every regulator with the adjusted meta-analysis.",
  S7 = "fgsea enrichment of positive- and negative-mode regulon targets on each cohort signature, and concordance with VIPER.")

wb <- createWorkbook()
bold <- createStyle(textDecoration = "bold")
hdr  <- createStyle(textDecoration = "bold", border = "bottom", fgFill = "#EDEDED")
sci  <- createStyle(numFmt = "0.00E+00")
dec  <- createStyle(numFmt = "0.000")
isp  <- function(n) grepl("(^|_|\\.)(p|pval|padj|P\\.Value|adj\\.P\\.Val|BH|p_value|p\\.value|Score|RRA_score|RRA_BH)($|_|\\.)|^p_|BH_|^p\\.", n)

# Índice
addWorksheet(wb, "Index")
idx <- data.frame(Table = paste0("Supplementary Table ", sapply(sheets, `[[`, 2)),
                  Sheet = sapply(sheets, `[[`, 1),
                  Rows = sapply(sheets, function(s) nrow(s[[4]])),
                  Description = sapply(sheets, `[[`, 3))
writeData(wb, "Index", "Supplementary Tables S1–S7 — Master regulator analysis of primary Sjögren's syndrome", startRow = 1)
addStyle(wb, "Index", bold, rows = 1, cols = 1)
writeData(wb, "Index", data.frame(Table = paste0("Supplementary Table ", names(legends)), Legend = legends),
          startRow = 3, headerStyle = hdr)
writeData(wb, "Index", idx, startRow = 3 + length(legends) + 3, headerStyle = hdr)
setColWidths(wb, "Index", cols = 1:4, widths = c(26, 20, 8, 120))

for (s in sheets) {
    nm <- s[[1]]; x <- s[[4]]
    addWorksheet(wb, nm)
    writeData(wb, nm, paste0("Supplementary Table ", s[[2]], ". ", s[[3]]), startRow = 1)
    addStyle(wb, nm, bold, rows = 1, cols = 1)
    # NA explícito: una celda vacía se omite del XML y algunos lectores recorren la fila
    writeData(wb, nm, x, startRow = 3, headerStyle = hdr, withFilter = TRUE,
              keepNA = TRUE, na.string = "NA")
    num <- which(sapply(x, function(v) is.numeric(v) && !all(v == round(v), na.rm = TRUE)))
    for (j in num) addStyle(wb, nm, if (isp(names(x)[j])) sci else dec,
                            rows = 4:(nrow(x) + 3), cols = j, gridExpand = TRUE)
    freezePane(wb, nm, firstActiveRow = 4, firstActiveCol = 2)
    setColWidths(wb, nm, cols = seq_along(x), widths = pmin(pmax(nchar(names(x)) + 2, 10), 40))
}
saveWorkbook(wb, out, overwrite = TRUE)
cat("hojas:", length(sheets) + 1, "\n")
