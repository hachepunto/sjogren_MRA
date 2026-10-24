suppressMessages({
    library(clusterProfiler)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Figura integradora de expresión diferencial (Figura 5). Solo dibuja: los
# resultados (incluido el ORA) vienen de compare_ED.R.
load("results/de/comparison_DE_results.rda")   # shared, sig_up, sig_dn, ego_de, ereact_de

BASE  <- 11   # tamaño base de letra
TITLE <- element_text(size = 14, hjust = 0.5)

# ── 1. Estado de cada gen compartido ──────────────────────────────────────────
# Replicados y discordantes: FDR < 0.05 en ambas cohortes (detalle en el pie de figura)
lv <- c("Replicated up", "Replicated down", "FDR < 0.05 one cohort", "Discordant",
        "Not significant")
shared$status <- lv[5]
shared$status[xor(shared$sig51, shared$sig84)] <- lv[3]
shared$status[shared$sig51 & shared$sig84 &
              sign(shared$logFC_51) != sign(shared$logFC_84)] <- lv[4]
shared$status[shared$concordant_up]   <- lv[1]
shared$status[shared$concordant_down] <- lv[2]
shared$status <- factor(shared$status, levels = lv)
cols_status <- setNames(c("#D73027", "#4575B4", "#FEE090", "#984EA3", "grey80"), lv)

# Etiquetas: replicados más alejados del origen en el plano de logFC
dist <- sqrt(shared$logFC_51^2 + shared$logFC_84^2)
top_far <- function(sel, n) shared[sel, ][head(order(-dist[sel]), n), ]
shared_sig <- rbind(top_far(shared$concordant_up, 20), top_far(shared$concordant_down, 8))

# ── 2. Panel A — scatter de logFC ─────────────────────────────────────────────
# Los no significativos se dibujan primero para que no tapen a los replicados
shared <- shared[order(shared$status, decreasing = TRUE), ]
pA <- ggplot(shared, aes(logFC_51, logFC_84, fill = status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(shape = 21, color = "white", stroke = 0.2, size = 1.6, alpha = 0.8) +
    geom_text_repel(data = shared_sig, aes(label = gene, color = status),
                    size = 3, max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, box.padding = 0.25, seed = 1,
                    show.legend = FALSE) +
    scale_fill_manual(values = cols_status, name = NULL) +
    scale_color_manual(values = cols_status, guide = "none") +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = expression(log[2]*FC ~ "—" ~ GSE51092 ~ "(n = 222)"),
         y = expression(log[2]*FC ~ "—" ~ GSE84844 ~ "(n = 60)"),
         title = "Cross-dataset differential expression") +
    theme_bw(base_size = BASE) +
    theme(plot.title = TITLE,
          legend.position = "inside",
          legend.position.inside = c(0.01, 0.99),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8.5))

# ── 3. Panel B — ORA (GO BP + Reactome) ──────────────────────────────────────
ora_go             <- as.data.frame(ego_de)[,    c("Description", "GeneRatio", "Count", "p.adjust")]
ora_go$Database    <- "GO BP"
ora_react          <- as.data.frame(ereact_de)[, c("Description", "GeneRatio", "Count", "p.adjust")]
ora_react$Database <- "Reactome"
ora_df             <- rbind(ora_go, ora_react)
ora_df$GeneRatioNum <- sapply(strsplit(ora_df$GeneRatio, "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))
ora_df$Description  <- factor(ora_df$Description,
                              levels = ora_df$Description[order(ora_df$GeneRatioNum)])

pB <- ggplot(ora_df, aes(x = GeneRatioNum, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatioNum, yend = Description, color = Database),
                 linewidth = 1.2) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = "white") +
    scale_color_manual(values = c("GO BP" = "#E08533", "Reactome" = "#7B3F9E"), name = NULL) +
    scale_fill_gradient(low = "#4575B4", high = "#FEE090", name = "FDR") +
    scale_size_continuous(name = "Count", range = c(2.5, 7)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    guides(color = guide_legend(order = 1),
           fill  = guide_colorbar(order = 2),
           size  = guide_legend(order = 3, override.aes = list(fill = "grey50"))) +
    labs(x = "Gene ratio", y = NULL,
         title = "Functional enrichment — replicated DE genes") +
    theme_bw(base_size = BASE) +
    theme(plot.title = TITLE,
          axis.text.y = element_text(size = 10),
          legend.position = "right",
          legend.box.spacing = unit(2, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = unit(4, "pt"))

# ── 4. Paneles C y D — top 30 replicados por |z Stouffer| ────────────────────
lollipop <- function(d, color, title) {
    d <- head(d[order(-abs(d$z_stouffer)), ], 30)
    d$logp <- -log10(d$p_stouffer)
    ggplot(d, aes(x = logp, y = reorder(gene, logp))) +
        geom_segment(aes(x = 0, xend = logp, yend = reorder(gene, logp)),
                     color = color, linewidth = 1) +
        geom_point(fill = color, color = "white", shape = 21, size = 2.6) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = expression(-log[10](italic(p)[Stouffer])), y = NULL, title = title) +
        theme_bw(base_size = BASE) +
        theme(plot.title = TITLE, axis.text.y = element_text(size = 9.5))
}
pC <- lollipop(shared[shared$concordant_up, ],   "#D73027", "Up-regulated in pSS\n(replicated)")
pD <- lollipop(shared[shared$concordant_down, ], "#4575B4", "Down-regulated in pSS\n(replicated)")

# ── 5. Composición: A y B en la columna ancha, C y D en la angosta; letras de
# panel independientes de los títulos
fig <- free(pA) + pB + pC + pD +   # free(): A no se alinea con las etiquetas largas de B
       plot_layout(design = "AC\nBD", widths = c(1, 0.4)) +
       plot_annotation(tag_levels = "A") &
       theme(plot.tag = element_text(size = 22, face = "bold"))

pdf("results/figures/DE_integrative_figure.pdf", width = 13, height = 12)
print(fig)
dev.off()
