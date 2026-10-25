# Carpetas de salida
invisible(lapply(c("results/integration"),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# Convergencia entre actividad (MRs de consenso, meta_MRs.R) y expresión
# (genes replicados, compare_ED.R / meta_DE.R): qué TFs son a la vez
# reguladores con z combinado concordante y genes replicados en la misma
# dirección.

mr <- read.delim("results/integration/MRs_meta_shared.txt")
de <- read.delim("results/integration/DE_meta_shared.txt")


m <- merge(mr, de[, c("gene", "direction", "replicated", "logFC_51", "adj.P.Val_51",
                      "logFC_84", "adj.P.Val_84", "z_stouffer", "BH_stouffer")],
           by.x = "Regulon", by.y = "gene", all.x = TRUE,
           suffixes = c("_MR", "_DE"))

m$DE_same_direction <- !is.na(m$replicated) & m$replicated &
    ((m$direction_MR == "activated" & m$direction_DE == "up") |
     (m$direction_MR == "repressed" & m$direction_DE == "down"))
m$measured_in_DE <- !is.na(m$direction_DE)

out <- m[order(m$p_stouffer), c("Regulon", "direction_MR", "z_stouffer_MR", "BH_stouffer_MR",
                                "consensus", "measured_in_DE", "direction_DE", "replicated",
                                "logFC_51", "adj.P.Val_51", "logFC_84", "adj.P.Val_84",
                                "z_stouffer_DE", "BH_stouffer_DE", "DE_same_direction")]
write.table(out, "results/integration/MR_DE_convergence.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Resumen: consenso (14) y 240 compartidos con signo concordante
summary_tab <- data.frame(
    set           = c("consensus_BH0.10", "shared_concordant_sign", "shared_all"),
    n_TFs         = c(sum(m$consensus), sum(m$direction_MR != "discordant"), nrow(m)),
    measured_in_DE = c(sum(m$consensus & m$measured_in_DE),
                       sum(m$direction_MR != "discordant" & m$measured_in_DE),
                       sum(m$measured_in_DE)),
    DE_replicated_same_dir = c(sum(m$consensus & m$DE_same_direction),
                               sum(m$direction_MR != "discordant" & m$DE_same_direction),
                               sum(m$DE_same_direction))
)
write.table(summary_tab, "results/integration/MR_DE_convergence_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
