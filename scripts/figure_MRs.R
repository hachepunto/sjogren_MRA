library(viper)

# Carpetas de salida
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Figura integradora de MRs. Requiere comparison_results.rda (compare_MRs.R) y
# MRs_meta_shared.txt (meta_MRs.R): correr ambos antes de este script.
load("results/integration/comparison_results.rda")   # shared

# ── 1. Consenso: Stouffer sobre NES con BH < 0.10 y mismo signo en ambas ────
meta <- read.delim("results/integration/MRs_meta_shared.txt")
meta$consensus <- meta$direction != "discordant" & meta$BH_stouffer < 0.10
sig_act <- meta$Regulon[meta$consensus & meta$direction == "activated"]
sig_rep <- meta$Regulon[meta$consensus & meta$direction == "repressed"]

# Colores en el scatter: rojo=activado, azul=reprimido, gris=resto
shared$color <- "grey80"
shared$color[shared$Regulon %in% sig_act] <- "#D73027"
shared$color[shared$Regulon %in% sig_rep] <- "#4575B4"
shared_sig <- shared[shared$Regulon %in% c(sig_act, sig_rep), ]

# Lollipop data: -log10(p Stouffer), consenso por dirección
lol_act <- meta[meta$Regulon %in% sig_act, ]
lol_act$logp <- -log10(lol_act$p_stouffer)
lol_act <- lol_act[order(lol_act$logp), ]

lol_rep <- meta[meta$Regulon %in% sig_rep, ]
lol_rep$logp <- -log10(lol_rep$p_stouffer)
lol_rep <- lol_rep[order(lol_rep$logp), ]

# Umbral de referencia: mayor p con BH < 0.10
p_thr <- max(meta$p_stouffer[meta$BH_stouffer < 0.10])

pdf("results/figures/MRs_integrative_figure.pdf", width = 14, height = 6)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.6, 1, 1))

# Panel A — Scatter NES
par(mar = c(5, 5, 4, 2))
plot(shared$NES_gse51092, shared$NES_84844,
     pch  = 21, bg = shared$color, col = "white", cex = 1.2,
     xlab = "NES — GSE51092 (ARACNE-AP, n=222)",
     ylab = "NES — GSE84844 (corto, n=60)",
     main = "A   Master Regulators cross-dataset")
abline(h = 0, v = 0, lty = 2, col = "grey60")
text(shared_sig$NES_gse51092, shared_sig$NES_84844,
     labels = shared_sig$Regulon,
     cex = 0.65, pos = 3, col = shared_sig$color)
legend("topleft",
       legend = c("Activated (consensus, FDR < 0.10)",
                  "Repressed (consensus, FDR < 0.10)", "Other"),
       pt.bg = c("#D73027", "#4575B4", "grey80"),
       pch = 21, col = "white", bty = "n", cex = 0.8)

# Función lollipop horizontal en base R
lollipop_h <- function(values, labels, color, main, xmax = NULL, vline) {
    n    <- length(values)
    xmax <- if (is.null(xmax)) max(values) * 1.1 else xmax
    par(mar = c(5, 6, 4, 2))
    plot(values, seq_len(n),
         xlim = c(0, xmax), ylim = c(0.5, n + 0.5),
         xlab = expression(-log[10](italic(p)[Stouffer])), ylab = "",
         main = main, yaxt = "n", pch = NA, bty = "l")
    axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.85)
    segments(x0 = 0, x1 = values, y0 = seq_len(n), col = color, lwd = 1.5)
    points(values, seq_len(n), pch = 21, bg = color, col = "white", cex = 1.6)
    abline(v = -log10(vline), lty = 2, col = "grey40")
}

xmax_shared <- max(c(lol_act$logp, lol_rep$logp)) * 1.1

# Panel B — consenso activados
lollipop_h(lol_act$logp, lol_act$Regulon, "#D73027",
           "B   Activated in pSS (consensus)", xmax_shared, p_thr)

# Panel C — consenso reprimidos
lollipop_h(lol_rep$logp, lol_rep$Regulon, "#4575B4",
           "C   Repressed in pSS (consensus)", xmax_shared, p_thr)

dev.off()
