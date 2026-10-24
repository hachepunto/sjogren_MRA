# Funciones compartidas por los scripts de MRA y expresión diferencial.

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

# Tabla de MRs significativos ordenada por |NES| (con leading edge si existe).
viper_mrsTopTable <- function(mrs, n = sum(mrs$es$p.value < p_threshold), p_threshold = 0.05) {
    sig <- mrs$es$p.value < p_threshold
    mrs_tops <- tibble::tibble(
        TF      = names(mrs$es$nes[sig]),
        nes     = mrs$es$nes[sig],
        p.value = mrs$es$p.value[sig]
    )
    mrs_tops <- mrs_tops[order(abs(mrs_tops$nes), decreasing = TRUE), ]
    mrs_tops <- mrs_tops[1:n, ]
    if (!is.null(mrs$ledge)) {
        mrs_tops$leading_edge <- sapply(mrs_tops$TF, function(tf) {
            le <- mrs$ledge[[tf]]
            if (is.null(le)) NA_character_ else paste(le, collapse = ",")
        })
    }
    return(mrs_tops)
}
