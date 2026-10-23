library(readr)
library(GSVA)
library(ggplot2)

# setwd("path/to/sjogren_MRA")  # set to repo root if not using an R project

# Carpetas de salida
invisible(lapply(c("results/validation", "results/figures"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ── 1. Processed gene-level expression matrix (GSE140161) ─────────────────────
# GSE140161_Processed_data.txt.gz: símbolos en la col 1 (sin nombre de columna)
# y 351 pacientes pSS en las cols 2-352. La cabecera tiene 351 nombres y cada
# fila 352 campos, así que read_tsv() corría las etiquetas una posición y
# perdía la última muestra; se lee la cabecera por separado. Sin controles.
hdr  <- strsplit(readLines(gzfile("data/GSE140161_Processed_data.txt.gz"), n = 1),
                 "\t")[[1]]
proc <- read.delim("data/GSE140161_Processed_data.txt.gz", header = FALSE,
                   skip = 1, check.names = FALSE)
stopifnot(ncol(proc) == length(hdr) + 1)
eset <- as.matrix(proc[, -1])
rownames(eset) <- proc[[1]]
colnames(eset) <- hdr
storage.mode(eset) <- "double"
stopifnot(!anyNA(eset))

# ── 2. Firma: genes up replicados (FDR < 0.05 en GSE51092 y GSE84844) ────────
sig_file <- "results/de/DE_replicated_up.txt"
suffix   <- "_replicated"
sig_up   <- read_tsv(sig_file, show_col_types = FALSE)
n_sig    <- nrow(sig_up)
gene_set <- sig_up$gene[sig_up$gene %in% rownames(eset)]

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
          paste0("results/validation/ssGSEA_GSE140161_scores", suffix, ".txt"))

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
                        length(gene_set), "/", n_sig, " signature genes present)")) +
    theme_bw(base_size = 11)

pdf(paste0("results/figures/ssGSEA_GSE140161", suffix, ".pdf"), width = 6, height = 5)
print(p)
dev.off()

# ── 5. Save ───────────────────────────────────────────────────────────────────
save(eset, scores, gene_set, n_sig, file = paste0("results/validation/ssGSEA_GSE140161", suffix, ".rda"))
