library(fgsea)
library(viper)

source("scripts/helpers.R")

# Enriquecimiento por rangos de los regulones completos sobre la firma de
# expresión de cada cohorte: complementa a VIPER con un
# estadístico independiente (fgsea multinivel) y p ajustada por BH sobre todos
# los regulones. Cada regulón se divide en blancos de modo positivo y negativo
# (fgsea no admite conjuntos con signo); el NES del conjunto positivo debe ir
# en el sentido de la actividad del TF y el del negativo en sentido contrario.

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
set.seed(1)
MINSIZE <- 10

split_regulon <- function(reg, genes) {
    sets <- list()
    for (tf in names(reg)) {
        tm <- reg[[tf]]$tfmode; tm <- tm[names(tm) %in% genes]
        pos <- names(tm)[tm > 0]; neg <- names(tm)[tm < 0]
        if (length(pos) >= MINSIZE) sets[[paste0(tf, "|pos")]] <- pos
        if (length(neg) >= MINSIZE) sets[[paste0(tf, "|neg")]] <- neg
    }
    sets
}
run_fgsea <- function(stats, sets, label) {
    stats <- stats[!is.na(stats)]
    res <- fgseaMultilevel(sets, stats, minSize = MINSIZE, maxSize = 2000,
                           eps = 0, nproc = 1)
    res <- as.data.frame(res[, c("pathway", "pval", "padj", "NES", "size")])
    res$Regulon <- sub("\\|.*", "", res$pathway)
    res$mode    <- sub(".*\\|", "", res$pathway)
    # NES orientado: el conjunto negativo se invierte para que ambos reflejen
    # la dirección de la actividad del TF
    res$NES_oriented <- ifelse(res$mode == "neg", -res$NES, res$NES)
    res$cohort <- label
    res[, c("cohort", "Regulon", "mode", "size", "NES", "NES_oriented", "pval", "padj")]
}

# ── 1. GSE51092: firma propia x regulones ARACNe-AP ──────────────────────────
e <- new.env(); load("results/mra/mra_GSE51092.rda", envir = e)
st51  <- signature_z(e$eset, e$disease)
sets51 <- split_regulon(e$regulon_gse51092, names(st51))
fg51  <- run_fgsea(st51, sets51, "GSE51092")
nes51 <- e$mrs_gse51092$es$nes
rm(e)

# ── 2. GSE84844: firma propia x regulones corto ──────────────────────────────
e <- new.env(); load("results/mra/mra_GSE84844.rda", envir = e)
nes84 <- e$mrs_gse84844$es$nes
st84 <- signature_z(e$eset, e$disease)
e2 <- new.env(); load("results/networks/regulon_GSE84844_corto.rda", envir = e2)
sets84 <- split_regulon(e2$regulon_gse84844, names(st84))
fg84  <- run_fgsea(st84, sets84, "GSE84844")
rm(e, e2)

fg <- rbind(fg51, fg84)
fg$NES_viper <- ifelse(fg$cohort == "GSE51092", nes51[fg$Regulon], nes84[fg$Regulon])

# ── 3. Concordancia fgsea vs VIPER (conjuntos positivos) ─────────────────────
ct <- function(d) {
    d <- d[d$mode == "pos" & !is.na(d$NES_viper), ]
    r <- cor.test(d$NES, d$NES_viper, method = "spearman", exact = FALSE)
    c(n = nrow(d), rho = unname(r$estimate), p = r$p.value)
}
concordance <- data.frame(cohort = c("GSE51092", "GSE84844"),
                          rbind(ct(fg51 |> transform(NES_viper = nes51[Regulon])),
                                ct(fg84 |> transform(NES_viper = nes84[Regulon]))))

# ── 4. Resumen por regulón: consenso y regulones de una sola red ─────────────
cons <- read_consensus()
wide <- reshape(fg[fg$mode == "pos", c("cohort", "Regulon", "NES", "padj")],
                idvar = "Regulon", timevar = "cohort", direction = "wide")
wide$consensus <- wide$Regulon %in% cons$Regulon
wide <- wide[order(!wide$consensus, wide$padj.GSE51092), ]

write.table(fg, "results/enrichment/fgsea_regulons_long.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wide, "results/enrichment/fgsea_regulons_wide.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concordance, "results/enrichment/fgsea_vs_viper_concordance.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(fg, wide, concordance, file = "results/enrichment/fgsea_regulons.rda")
