library(viper)

source("scripts/helpers.R")

# Validación con interactoma externo curado: VIPER sobre las firmas pSS vs
# control de ambas cohortes usando los regulones de CollecTRI (OmniPath;
# resources/collectri_regulons_human.tsv, descargado de
# https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=1,
# complejos excluidos, mor = -1 si solo inhibición, +1 en otro caso).
# Los regulones curados son más pequeños que los inferidos, por lo que se usa
# minsize = 10 (en lugar de 25) y se reporta el tamaño efectivo por TF.

dir.create("results/robustness", recursive = TRUE, showWarnings = FALSE)
set.seed(1)
MINSIZE <- 10

# ── 1. Regulón CollecTRI en formato viper ────────────────────────────────────
ct <- read.delim("resources/collectri_regulons_human.tsv")
build_regulon <- function(ct, genes) {
    ct <- ct[ct$target %in% genes & ct$source != ct$target, ]
    sp <- split(ct, ct$source)
    reg <- lapply(sp, function(d) {
        d <- d[!duplicated(d$target), ]
        list(tfmode = setNames(d$mor, d$target),
             likelihood = rep(1, nrow(d)))
    })
    reg <- reg[sapply(reg, function(r) length(r$tfmode)) >= MINSIZE]
    class(reg) <- "regulon"
    reg
}

# ── 2. Firmas y VIPER por cohorte ────────────────────────────────────────────

e <- new.env(); load("results/mra/mra_GSE51092.rda", envir = e)
reg51 <- build_regulon(ct, rownames(e$eset))
ct51  <- run_msviper(e$eset, e$disease, reg51, minsize = MINSIZE)
size51 <- sapply(reg51, function(r) length(r$tfmode))

e <- new.env(); load("results/mra/mra_GSE84844.rda", envir = e)
reg84 <- build_regulon(ct, rownames(e$eset))
ct84  <- run_msviper(e$eset, e$disease, reg84, minsize = MINSIZE)
size84 <- sapply(reg84, function(r) length(r$tfmode))
rm(e)

# ── 3. Tabla, meta-análisis Stouffer y comparación con el consenso ───────────
t51 <- data.frame(Regulon = names(ct51$es$nes), size_51092 = size51[names(ct51$es$nes)],
                  NES_ct_51092 = ct51$es$nes, p_ct_51092 = ct51$es$p.value)
t84 <- data.frame(Regulon = names(ct84$es$nes), size_84844 = size84[names(ct84$es$nes)],
                  NES_ct_84844 = ct84$es$nes, p_ct_84844 = ct84$es$p.value)
res <- merge(t51, t84, by = "Regulon", all = FALSE)
res$same_sign  <- sign(res$NES_ct_51092) == sign(res$NES_ct_84844)
res$z_stouffer <- (res$NES_ct_51092 + res$NES_ct_84844) / sqrt(2)
res$p_stouffer <- 2 * pnorm(-abs(res$z_stouffer))
res$BH_stouffer <- p.adjust(res$p_stouffer, "BH")

# Consenso de las redes inferidas (meta_MRs.R) y cobertura en CollecTRI
cons <- read.delim("results/integration/MRs_meta_shared.txt")
cons$in_collectri <- cons$Regulon %in% ct$source
cons$collectri_size_51092 <- size51[cons$Regulon]
cons$collectri_size_84844 <- size84[cons$Regulon]
cons$NES_ct_51092 <- res$NES_ct_51092[match(cons$Regulon, res$Regulon)]
cons$NES_ct_84844 <- res$NES_ct_84844[match(cons$Regulon, res$Regulon)]
cons$z_ct_stouffer <- res$z_stouffer[match(cons$Regulon, res$Regulon)]
cons$BH_ct_stouffer <- res$BH_stouffer[match(cons$Regulon, res$Regulon)]
res$consensus_inferred <- res$Regulon %in% cons$Regulon[cons$consensus]

# Concordancia global entre NES inferido y NES CollecTRI (TFs en ambos)
both <- merge(cons[, c("Regulon", "z_stouffer")], res[, c("Regulon", "z_stouffer")],
              by = "Regulon", suffixes = c("_inferred", "_collectri"))
rho <- cor.test(both$z_stouffer_inferred, both$z_stouffer_collectri,
                method = "spearman", exact = FALSE)
concordance <- data.frame(n_TFs = nrow(both), rho = unname(rho$estimate),
                          p = rho$p.value)

res <- res[order(res$p_stouffer), ]
write.table(res,  "results/robustness/viper_collectri.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cons, "results/robustness/consensus_vs_collectri.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(concordance, "results/robustness/viper_collectri_concordance.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
save(ct51, ct84, res, cons, concordance, file = "results/robustness/viper_collectri.rda")
