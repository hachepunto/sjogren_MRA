# Cifras de composición celular que se citan en la Sección 2.3 de los
# Resultados Suplementarios, tomadas tal cual de deconv_composition_tests.txt
# (BH dentro de cada cohorte sobre los 16 tipos celulares de ambos métodos).
comp <- read.delim("results/robustness/deconv_composition_tests.txt", check.names = FALSE)

cited <- data.frame(
    method    = c("xCell", "xCell", "MCPcounter", "MCPcounter", "xCell", "xCell",
                  "MCPcounter", "MCPcounter"),
    cell_type = c("CD8+ T-cells", "B-cells", "B lineage", "Monocytic lineage", "NK cells",
                  "Neutrophils", "NK cells", "Neutrophils"),
    statement = c(rep("cited", 4), rep("no difference", 4))
)
out <- merge(cited, comp, by = c("method", "cell_type"))
stopifnot(nrow(out) == 2 * nrow(cited))
out <- out[order(out$statement, out$method, out$cell_type, out$cohort),
           c("cell_type", "cohort", "method", "median_pSS", "median_control", "p", "BH", "statement")]
print(out, row.names = FALSE)
write.table(out, "results/robustness/composition_reported.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
