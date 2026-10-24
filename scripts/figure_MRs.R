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

# Lollipop: -log10(p Stouffer) de todos los MRs de consenso, color por dirección
lol <- meta[meta$Regulon %in% c(sig_act, sig_rep), ]
lol$logp  <- -log10(lol$p_stouffer)
lol$color <- ifelse(lol$direction == "activated", "#D73027", "#4575B4")
lol <- lol[order(lol$logp), ]

# Umbral de referencia: mayor p con BH < 0.10
p_thr <- max(meta$p_stouffer[meta$BH_stouffer < 0.10])

pdf("results/figures/MRs_integrative_figure.pdf", width = 11, height = 6)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.6, 1))

# Panel A — Scatter NES
par(mar = c(5, 5, 4, 2))
plot(shared$NES_gse51092, shared$NES_84844,
     pch  = 21, bg = shared$color, col = "white", cex = 1.2,
     ylim = range(shared$NES_84844) + c(0, 0.3),   # espacio para las etiquetas de arriba
     xlab = "NES — GSE51092 (ARACNE-AP, n=222)",
     ylab = "NES — GSE84844 (corto, n=60)",
     main = "A   Master Regulators cross-dataset")
abline(h = 0, v = 0, lty = 2, col = "grey60")
# Posición de cada etiqueta (1 abajo, 2 izq., 3 arriba, 4 der.) para evitar
# que se encimen en el cuadrante superior derecho
lab_pos <- c(EEA1 = 2, TFEC = 3, ZNF184 = 4, ZNF140 = 4, ZNF227 = 2, ZHX1 = 1,
             ZNF654 = 3, ELF1 = 3, NFE2L3 = 3, PLSCR1 = 1)
lab_dy  <- c(ZNF140 = 0.08, ZNF654 = -0.04)   # ajuste vertical fino
pos <- ifelse(shared_sig$Regulon %in% names(lab_pos), lab_pos[shared_sig$Regulon], 3)
dy  <- ifelse(shared_sig$Regulon %in% names(lab_dy), lab_dy[shared_sig$Regulon], 0)
text(shared_sig$NES_gse51092, shared_sig$NES_84844 + dy,
     labels = shared_sig$Regulon,
     cex = 0.65, pos = pos, offset = 0.4, col = shared_sig$color)
legend("bottomright",
       legend = c("Activated (consensus, FDR < 0.10)",
                  "Repressed (consensus, FDR < 0.10)", "Other"),
       pt.bg = c("#D73027", "#4575B4", "grey80"),
       pch = 21, col = "white", bty = "n", cex = 0.8)

# Panel B — lollipop horizontal de los MRs de consenso
n <- nrow(lol)
par(mar = c(5, 6, 4, 2))
plot(lol$logp, seq_len(n),
     xlim = c(0, max(lol$logp) * 1.1), ylim = c(0.5, n + 3),
     xlab = expression(-log[10](italic(p)[Stouffer])), ylab = "",
     main = "B   Consensus MRs", yaxt = "n", pch = NA, bty = "l")
axis(2, at = seq_len(n), labels = lol$Regulon, las = 1, cex.axis = 0.85)
segments(x0 = 0, x1 = lol$logp, y0 = seq_len(n), col = lol$color, lwd = 1.5)
points(lol$logp, seq_len(n), pch = 21, bg = lol$color, col = "white", cex = 1.6)
abline(v = -log10(p_thr), lty = 2, col = "grey40")
legend("topleft", legend = c("Activated in pSS", "Repressed in pSS"),
       pt.bg = c("#D73027", "#4575B4"), pch = 21, col = "white",
       pt.cex = 1.4, bty = "n", cex = 0.8)

dev.off()
