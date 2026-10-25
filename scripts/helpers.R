# Funciones compartidas por los scripts de MRA, expresión diferencial y robustez.

# Colapsa sondas a genes: por símbolo conserva la sonda de mayor mediana.
collapse_probes <- function(eset, symbols) {
    keep    <- !is.na(symbols) & symbols != "" & !grepl("///", symbols)
    eset    <- eset[keep, ]
    symbols <- symbols[keep]
    eset <- do.call(rbind, lapply(
        split(seq_len(nrow(eset)), symbols),
        function(idx) {
            if (length(idx) == 1) return(eset[idx, , drop = FALSE])
            eset[idx[which.max(matrixStats::rowMedians(eset[idx, , drop = FALSE]))], , drop = FALSE]
        }
    ))
    rownames(eset) <- sort(unique(symbols))
    eset
}

# Firma pSS vs control: t de Welch por gen convertido a z con su signo.
signature_z <- function(eset, disease) {
    cases <- which(disease == "pSS"); ctrls <- which(disease == "control")
    s <- viper::rowTtest(eset[, cases], eset[, ctrls])
    (qnorm(s$p.value / 2, lower.tail = FALSE) * sign(s$statistic))[, 1]
}

# msVIPER con la firma pSS vs control y un nulo de 1000 permutaciones de muestras.
run_msviper <- function(eset, disease, regulon, minsize = 25, per = 1000) {
    cases <- which(disease == "pSS"); ctrls <- which(disease == "control")
    nullm <- viper::ttestNull(eset[, cases], eset[, ctrls], per = per, repos = TRUE,
                              verbose = FALSE)
    viper::msviper(signature_z(eset, disease), regulon, nullm, minsize = minsize,
                   verbose = FALSE)
}

# Símbolos HGNC → Entrez (los símbolos sin correspondencia se descartan).
sym2entrez <- function(genes) {
    clusterProfiler::bitr(as.character(genes), fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE)$ENTREZID
}

# Reguladores de consenso (columna escrita una sola vez por meta_MRs.R).
read_consensus <- function() {
    cons <- read.delim("results/integration/MRs_meta_shared.txt")
    cons[cons$consensus, ]
}

# Tabla de MRs significativos ordenada por |NES| (con leading edge si existe).
viper_mrsTopTable <- function(mrs, n = sum(mrs$es$p.value < p_threshold), p_threshold = 0.05) {
    sig <- mrs$es$p.value < p_threshold
    mrs_tops <- tibble::tibble(
        TF      = names(mrs$es$nes[sig]),
        nes     = mrs$es$nes[sig],
        p.value = mrs$es$p.value[sig]
    )
    mrs_tops <- mrs_tops[order(abs(mrs_tops$nes), decreasing = TRUE), ]
    mrs_tops <- mrs_tops[seq_len(min(n, nrow(mrs_tops))), ]   # sin filas NA si hay menos de n
    if (!is.null(mrs$ledge)) {
        mrs_tops$leading_edge <- sapply(mrs_tops$TF, function(tf) {
            le <- mrs$ledge[[tf]]
            if (is.null(le)) NA_character_ else paste(le, collapse = ",")
        })
    }
    return(mrs_tops)
}
