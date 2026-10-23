# Reconstrucción de la matriz de expresión de GSE84844 con el mismo
# preprocesamiento que mra_GSE84844.R (colapso de sondas por mediana máxima).
# Se cachea en results/mra/eset_GSE84844.rda para los análisis de robustez.
get_eset_gse84844 <- function(cache = "results/mra/eset_GSE84844.rda") {
    if (file.exists(cache)) { e <- new.env(); load(cache, envir = e); return(e) }
    suppressMessages({ library(GEOquery); library(matrixStats) })
    gse  <- getGEO(filename = "data/GSE84844_series_matrix.txt",
                   GSEMatrix = TRUE, AnnotGPL = TRUE)
    eset <- exprs(gse); sym <- fData(gse)[["Gene symbol"]]
    keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
    eset <- eset[keep, ]; sym <- sym[keep]
    eset <- do.call(rbind, lapply(split(seq_len(nrow(eset)), sym), function(idx) {
        if (length(idx) == 1) return(eset[idx, , drop = FALSE])
        eset[idx[which.max(rowMedians(eset[idx, , drop = FALSE]))], , drop = FALSE]
    }))
    rownames(eset) <- sort(unique(sym))
    disease <- ifelse(pData(gse)[["disease:ch1"]] == "Healthy control", "control", "pSS")
    stopifnot(table(disease)[["control"]] == 30, table(disease)[["pSS"]] == 30)
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    save(eset, disease, file = cache)
    e <- new.env(); e$eset <- eset; e$disease <- disease; e
}
