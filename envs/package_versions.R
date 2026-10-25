# Escribe package_versions.tsv con las versiones del ambiente de R activo: los
# paquetes que cargan los scripts y las bases de anotación que deciden el ORA.
# Uso: Rscript envs/package_versions.R
pkgs <- c("viper", "corto", "limma", "GEOquery", "matrixStats", "RobustRankAggreg",
          "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "GO.db", "reactome.db",
          "enrichplot", "GSVA", "GSEABase", "fgsea", "xCell", "MCPcounter",
          "ggplot2", "ggrepel", "patchwork", "openxlsx", "readr", "tibble")
src <- function(p) {
    d <- packageDescription(p)
    if (!is.null(d$RemoteSha)) sprintf("GitHub %s/%s@%s", d$RemoteUsername, d$RemoteRepo, d$RemoteSha)
    else if (!is.null(d$biocViews)) "Bioconductor" else "CRAN"
}
tab <- data.frame(Package = c("R", "Bioconductor", pkgs),
                  Version = c(paste(R.version$major, R.version$minor, sep = "."),
                              as.character(BiocManager::version()),
                              sapply(pkgs, function(p) as.character(packageVersion(p)))),
                  Source  = c("", "", sapply(pkgs, src)))
write.table(tab, "package_versions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, row.names = FALSE)
