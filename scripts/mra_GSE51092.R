library(limma)
library(viper)
library(tibble)
library(readr)

source("scripts/helpers.R")

invisible(lapply(c("results/mra", "results/de", "results/figures"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

load("results/mra/eset_GSE51092.rda")   # eset, disease (prep_GSE51092.R)

# ── 1. Regulon desde ARACNe-AP ────────────────────────────────────────────────
regulon_gse51092 <- aracne2regulon("results/networks/sjs_network.txt", eset)
cat("Regulon:", length(regulon_gse51092), "TFs\n")

# ── 2. Firma molecular ────────────────────────────────────────────────────────
Ncases    <- which(disease == "pSS")
Ncontrols <- which(disease == "control")

signature <- signature_z(eset, disease)

# ── 3. Null model ─────────────────────────────────────────────────────────────
set.seed(1)   # permutaciones de ttestNull
nullmodel <- ttestNull(eset[, Ncases], eset[, Ncontrols],
                       per = 1000, repos = TRUE)

# ── 4. Master Regulator Analysis ──────────────────────────────────────────────
mrs_gse51092 <- msviper(signature, regulon_gse51092, nullmodel)
print(summary(mrs_gse51092))
mrs_gse51092 <- ledge(mrs_gse51092)

# ── 5. Resultados MRA ─────────────────────────────────────────────────────────
top_mrs <- viper_mrsTopTable(mrs_gse51092, n = 100)
write_tsv(top_mrs, "results/mra/top100_MRs_GSE51092.txt")

pdf("results/figures/Top10_MRs_GSE51092.pdf", width = 6, height = 7)
plot(mrs_gse51092, mrs = 10, cex = 0.7)
dev.off()

# ── 6. Expresión diferencial (limma) ──────────────────────────────────────────
disease_f      <- factor(disease, levels = c("control", "pSS"))
design         <- model.matrix(~ disease_f)
fit            <- lmFit(eset, design)
fit            <- eBayes(fit)
de_gse51092   <- topTable(fit, coef = 2, number = Inf,
                           adjust.method = "BH", sort.by = "P")
de_gse51092   <- tibble::rownames_to_column(de_gse51092, var = "gene")
write_tsv(de_gse51092, "results/de/DE_limma_GSE51092.txt")
cat("DE: ", sum(de_gse51092$adj.P.Val < 0.05), "genes significativos (FDR < 0.05)\n")

save(mrs_gse51092, top_mrs, regulon_gse51092, eset, disease, de_gse51092,
     file = "results/mra/mra_GSE51092.rda")
