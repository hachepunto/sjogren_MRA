library(readr)
library(GSVA)
library(ggplot2)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# ── 1. Processed gene-level expression matrix (GSE140161) ─────────────────────
# GSE140161_Processed_data.txt.gz: gene symbols in col 1 (no column header),
# 350 pSS patient samples in cols 2-351. No healthy controls in this dataset.
proc    <- read_tsv("data/GSE140161_Processed_data.txt.gz",
                    show_col_types = FALSE)
eset    <- as.matrix(proc[, -1])
rownames(eset) <- proc[[1]]
storage.mode(eset) <- "double"

# ── 2. IFN-up gene set (RRA 636-gene upregulated signature) ───────────────────
rra_up   <- read_tsv("DE_RRA_upregulated.txt", show_col_types = FALSE)
gene_set <- rra_up$Name[rra_up$Name %in% rownames(eset)]  # 576 / 636 present

# ── 3. ssGSEA (Barbie et al. 2009, per-sample enrichment) ────────────────────
param  <- ssgseaParam(exprData  = eset,
                      geneSets  = list(IFN_up = gene_set),
                      alpha     = 0.25,
                      normalize = TRUE)
es     <- gsva(param, verbose = FALSE)
scores <- as.numeric(es["IFN_up", ])
names(scores) <- colnames(eset)

write_tsv(data.frame(sample      = names(scores),
                     ssGSEA_score = scores),
          "ssGSEA_GSE140161_scores.txt")

print(summary(scores))

# ── 4. Density figure ─────────────────────────────────────────────────────────
df_plot <- data.frame(score = scores)

p <- ggplot(df_plot, aes(x = score)) +
    geom_density(fill = "#4393C3", alpha = 0.65, linewidth = 0.7) +
    geom_vline(xintercept = median(scores), linetype = "dashed",
               color = "gray30", linewidth = 0.6) +
    annotate("text", x = median(scores), y = Inf,
             label = paste0("Median = ", round(median(scores), 3)),
             vjust = 1.5, hjust = -0.05, size = 3.5, color = "gray30") +
    labs(x = "ssGSEA enrichment score — IFN-up signature",
         y = "Density",
         title = paste0("IFN-up signature enrichment in GSE140161\n",
                        "(n = ", ncol(eset), " pSS patients, ",
                        length(gene_set), "/636 signature genes present)")) +
    theme_bw(base_size = 11)

pdf("ssGSEA_GSE140161.pdf", width = 6, height = 5)
print(p)
dev.off()

# ── 5. Save ───────────────────────────────────────────────────────────────────
save(eset, scores, gene_set, file = "ssGSEA_GSE140161.rda")
