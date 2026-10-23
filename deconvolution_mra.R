library(xCell)
library(MCPcounter)
library(limma)
library(viper)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# Análisis de sensibilidad a la composición celular (Review 4, punto M1).
# Sangre total: la firma pSS vs control mezcla cambios por célula y cambios de
# proporciones. Se estiman proporciones relativas con xCell (y MCP-counter como
# segundo estimador), se contrasta pSS vs control por tipo celular, y se repite
# la MRA sobre la expresión ajustada por las puntuaciones xCell de seis
# poblaciones. Advertencia: el IFN tipo I contribuye a la linfopenia en pSS,
# así que el ajuste puede eliminar biología real; se reportan ambas versiones.

dir.create("results/robustness", recursive = TRUE, showWarnings = FALSE)
set.seed(1)

CELLS <- c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "NK cells",
           "Monocytes", "Neutrophils")

deconv <- function(eset) {
    xc <- xCellAnalysis(eset, rnaseq = FALSE, parallel.sz = 1)
    mc <- MCPcounter.estimate(eset, featuresType = "HUGO_symbols")
    list(xcell = xc, mcp = mc)
}
test_composition <- function(scores, disease, label, method) {
    do.call(rbind, lapply(rownames(scores), function(ct) {
        x <- scores[ct, disease == "pSS"]; y <- scores[ct, disease == "control"]
        w <- wilcox.test(x, y)
        data.frame(cohort = label, method = method, cell_type = ct,
                   median_pSS = median(x), median_control = median(y),
                   p = w$p.value)
    }))
}
adjusted_mra <- function(eset, disease, regulon, covars, per = 1000) {
    covars <- t(scale(t(covars)))                  # tipos celulares x muestras
    adj    <- removeBatchEffect(eset, covariates = t(covars))
    cases  <- which(disease == "pSS"); ctrls <- which(disease == "control")
    sig    <- rowTtest(adj[, cases], adj[, ctrls])
    sig    <- (qnorm(sig$p.value / 2, lower.tail = FALSE) * sign(sig$statistic))[, 1]
    nullm  <- ttestNull(adj[, cases], adj[, ctrls], per = per, repos = TRUE,
                        verbose = FALSE)
    msviper(sig, regulon, nullm, minsize = 25, verbose = FALSE)
}

run_cohort <- function(eset, disease, regulon, mrs_own, label) {
    dc  <- deconv(eset)
    comp <- rbind(test_composition(dc$xcell[CELLS, ], disease, label, "xCell"),
                  test_composition(dc$mcp, disease, label, "MCPcounter"))
    comp$BH <- p.adjust(comp$p, "BH")
    adj <- adjusted_mra(eset, disease, regulon, dc$xcell[CELLS, ])
    tab <- data.frame(Regulon = names(mrs_own$es$nes),
                      NES_unadj = mrs_own$es$nes, p_unadj = mrs_own$es$p.value)
    tab$NES_adj <- adj$es$nes[tab$Regulon]
    tab$p_adj   <- adj$es$p.value[tab$Regulon]
    write.table(t(dc$xcell), sprintf("results/robustness/xcell_%s.txt", label),
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(t(dc$mcp), sprintf("results/robustness/mcpcounter_%s.txt", label),
                sep = "\t", quote = FALSE, col.names = NA)
    list(comp = comp, tab = tab, mrs_adj = adj)
}

# ── 1. GSE51092 ───────────────────────────────────────────────────────────────
e <- new.env(); load("results/mra/mra_GSE51092.rda", envir = e)
r51 <- run_cohort(e$eset, e$disease, e$regulon_gse51092, e$mrs_gse51092, "GSE51092")
rm(e)

# ── 2. GSE84844 ───────────────────────────────────────────────────────────────
e <- new.env(); load("results/mra/mra_GSE84844.rda", envir = e)
mrs84 <- e$mrs_gse84844
if (is.null(e$eset)) { source("helpers_eset.R"); e <- get_eset_gse84844() }
e2 <- new.env(); load("results/networks/regulon_GSE84844_corto.rda", envir = e2)
r84 <- run_cohort(e$eset, e$disease, e2$regulon_gse84844, mrs84, "GSE84844")
rm(e, e2)

# ── 3. Resumen: concordancia ajustado vs no ajustado y consenso ajustado ─────
ct <- function(t) {
    ok <- complete.cases(t$NES_unadj, t$NES_adj)
    r  <- cor.test(t$NES_unadj[ok], t$NES_adj[ok], method = "spearman", exact = FALSE)
    c(n = sum(ok), rho = unname(r$estimate), p = r$p.value)
}
concordance <- data.frame(cohort = c("GSE51092", "GSE84844"),
                          rbind(ct(r51$tab), ct(r84$tab)))

tab <- merge(r51$tab, r84$tab, by = "Regulon", suffixes = c("_51092", "_84844"))
stouffer <- function(a, b) {
    z <- (a + b) / sqrt(2); p <- 2 * pnorm(-abs(z))
    data.frame(z = z, p = p, BH = p.adjust(p, "BH"), same_sign = sign(a) == sign(b))
}
m_un <- stouffer(tab$NES_unadj_51092, tab$NES_unadj_84844); names(m_un) <- paste0("unadj_", names(m_un))
m_ad <- stouffer(tab$NES_adj_51092,   tab$NES_adj_84844);   names(m_ad) <- paste0("adj_",   names(m_ad))
tab <- cbind(tab, m_un, m_ad)
cons <- read.delim("results/integration/MRs_meta_shared.txt")
tab$consensus <- tab$Regulon %in% cons$Regulon[cons$BH_stouffer < 0.10 & cons$direction != "discordant"]
tab <- tab[order(tab$adj_p), ]

write.table(rbind(r51$comp, r84$comp), "results/robustness/deconv_composition_tests.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab, "results/robustness/mra_adjusted_composition.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concordance, "results/robustness/mra_adjusted_concordance.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(r51, r84, tab, concordance, file = "results/robustness/deconvolution_mra.rda")
