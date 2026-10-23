library(viper)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# Validación cruzada red/firma (Reviewer 1, punto 4): la actividad de cada TF
# se estima con la red inferida en la OTRA cohorte, de modo que la red y la
# firma de expresión nunca provienen de las mismas muestras.
#   - firma GSE84844 (pSS vs control) x red ARACNE-AP de GSE51092
#   - firma GSE51092 (pSS vs control) x red corto de GSE84844
# Se compara con la actividad obtenida con la red propia (mra_*.rda).

dir.create("results/robustness", recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# ── 1. Cargar esets, firmas y redes ──────────────────────────────────────────
e <- new.env()
load("results/mra/mra_GSE51092.rda", envir = e)   # mrs_gse51092, regulon_gse51092, eset, disease
eset51 <- e$eset; dis51 <- e$disease; reg51 <- e$regulon_gse51092; mrs51 <- e$mrs_gse51092
e <- new.env()
load("results/mra/mra_GSE84844.rda", envir = e)   # mrs_gse84844 (+ eset, disease si es reciente)
mrs84 <- e$mrs_gse84844
if (is.null(e$eset)) { source("helpers_eset.R"); e <- get_eset_gse84844() }
eset84 <- e$eset; dis84 <- e$disease
e <- new.env()
load("results/networks/regulon_GSE84844_corto.rda", envir = e)   # regulon_gse84844
reg84 <- e$regulon_gse84844
rm(e)

run_msviper <- function(eset, disease, regulon, per = 1000) {
    cases <- which(disease == "pSS"); ctrls <- which(disease == "control")
    sig   <- rowTtest(eset[, cases], eset[, ctrls])
    sig   <- (qnorm(sig$p.value / 2, lower.tail = FALSE) * sign(sig$statistic))[, 1]
    nullm <- ttestNull(eset[, cases], eset[, ctrls], per = per, repos = TRUE,
                       verbose = FALSE)
    msviper(sig, regulon, nullm, minsize = 25, verbose = FALSE)
}

# ── 2. VIPER cruzado ──────────────────────────────────────────────────────────
cross84 <- run_msviper(eset84, dis84, reg51)   # firma 84844 x red 51092
cross51 <- run_msviper(eset51, dis51, reg84)   # firma 51092 x red 84844

tab <- function(mrs, label) {
    data.frame(Regulon = names(mrs$es$nes),
               nes = mrs$es$nes, p = mrs$es$p.value,
               stringsAsFactors = FALSE) |>
        setNames(c("Regulon", paste0("NES_", label), paste0("p_", label)))
}
own51   <- tab(mrs51,   "own_51092")      # firma 51092 x red 51092
own84   <- tab(mrs84,   "own_84844")      # firma 84844 x red 84844
x84     <- tab(cross84, "cross_84844")    # firma 84844 x red 51092
x51     <- tab(cross51, "cross_51092")    # firma 51092 x red 84844

res <- Reduce(function(a, b) merge(a, b, by = "Regulon", all = TRUE),
              list(own51, own84, x84, x51))

# ── 3. Concordancia propia vs cruzada ────────────────────────────────────────
# (a) misma firma, red distinta: ¿la red cambia la estimación de actividad?
# (b) misma red, firma distinta: replicación de la actividad con red fija.
ct <- function(x, y) {
    ok <- complete.cases(x, y)
    r  <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    c(n = sum(ok), rho = unname(r$estimate), p = r$p.value)
}
concordance <- rbind(
    "sig84844: red84844 vs red51092"   = ct(res$NES_own_84844, res$NES_cross_84844),
    "sig51092: red51092 vs red84844"   = ct(res$NES_own_51092, res$NES_cross_51092),
    "red51092: sig51092 vs sig84844"   = ct(res$NES_own_51092, res$NES_cross_84844),
    "red84844: sig84844 vs sig51092"   = ct(res$NES_own_84844, res$NES_cross_51092)
)
concordance <- data.frame(comparison = rownames(concordance), concordance,
                          row.names = NULL)

# ── 4. Consenso con red fija (sin dependencia red/firma) ─────────────────────
# Red ARACNE-AP (GSE51092) aplicada a ambas firmas: Stouffer de NES propio
# (51092) y cruzado (84844); idem para la red corto.
meta <- function(z1, z2) {
    z <- (z1 + z2) / sqrt(2); p <- 2 * pnorm(-abs(z))
    data.frame(z = z, p = p, BH = p.adjust(p, "BH"),
               same_sign = sign(z1) == sign(z2))
}
m_aracne <- meta(res$NES_own_51092, res$NES_cross_84844)
m_corto  <- meta(res$NES_cross_51092, res$NES_own_84844)
names(m_aracne) <- paste0("aracne_", names(m_aracne))
names(m_corto)  <- paste0("corto_",  names(m_corto))
res <- cbind(res, m_aracne, m_corto)

# ── 5. Recuperación de los MRs de consenso ───────────────────────────────────
cons <- read.delim("results/integration/MRs_meta_shared.txt")
cons <- cons[cons$BH_stouffer < 0.10, c("Regulon", "direction", "z_stouffer")]
res$consensus <- res$Regulon %in% cons$Regulon
res <- res[order(res$aracne_p), ]

write.table(res, "results/robustness/viper_cross_network.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concordance, "results/robustness/viper_cross_network_concordance.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(cross84, cross51, res, concordance,
     file = "results/robustness/viper_cross_network.rda")
