# Tabla 1 del texto principal: resumen de los 14 reguladores de
# consenso que reúne actividad (meta_MRs.R), expresión del propio transcrito
# (meta_DE.R), ORA del regulón (ora_regulons.R) y cobertura en CollecTRI
# (viper_collectri.R). Una fila por regulador, ordenada por |z| descendente.
dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

load("results/mra/mra_GSE51092.rda")                 # regulon_gse51092
mr <- read.delim("results/integration/MRs_meta_shared.txt")
de <- read.delim("results/integration/DE_meta_shared.txt")
ct <- read.delim("results/robustness/consensus_vs_collectri.txt")
go <- read.delim("results/enrichment/ORA_GO_BP.txt")
re <- read.delim("results/enrichment/ORA_Reactome.txt")

cons <- mr[mr$consensus, ]
cons <- cons[order(-abs(cons$z_stouffer)), ]
tf   <- cons$Regulon

# Término más significativo por regulador; empates por p y luego alfabético
top_term <- function(ora) {
    ora <- ora[order(ora$Cluster, ora$p.adjust, ora$pvalue, ora$Description), ]
    first <- ora[!duplicated(ora$Cluster), ]
    setNames(first$Description, first$Cluster)[tf]
}
top_go <- top_term(go)
top_re <- top_term(re)
both   <- rbind(transform(go, db = "GO BP"), transform(re, db = "Reactome"))
both   <- both[order(both$Cluster, both$p.adjust, both$pvalue, both$Description), ]
best   <- both[!duplicated(both$Cluster), ]
best   <- best[match(tf, best$Cluster), ]

dx <- de[match(tf, de$gene), ]
cx <- ct[match(tf, ct$Regulon), ]
# Blancos curados medidos por cohorte; viper_collectri.R solo guarda regulones con >= 10
fmt_ct  <- function(n) ifelse(is.na(n), "<10", n)
curated <- ifelse(cx$in_collectri,
                  paste(fmt_ct(cx$collectri_size_51092), fmt_ct(cx$collectri_size_84844), sep = " / "),
                  "absent")

tab <- data.frame(
    regulator             = tf,
    direction             = cons$direction,
    regulon_size_ARACNe   = sapply(regulon_gse51092[tf], function(r) length(r$tfmode)),
    NES_GSE51092          = cons$NES_gse51092,
    NES_GSE84844          = cons$NES_84844,
    z_stouffer            = cons$z_stouffer,
    BH_stouffer           = cons$BH_stouffer,
    logFC_51              = dx$logFC_51,
    adjP_51               = dx$adj.P.Val_51,
    logFC_84              = dx$logFC_84,
    adjP_84               = dx$adj.P.Val_84,
    transcript_replicated = as.integer(!is.na(dx$replicated) & dx$replicated),
    top_term              = ifelse(is.na(best$Description), "none", best$Description),
    top_term_db           = ifelse(is.na(best$db), "", best$db),
    top_term_padj         = best$p.adjust,
    top_GO_BP             = ifelse(is.na(top_go), "none", top_go),
    top_Reactome          = ifelse(is.na(top_re), "none", top_re),
    collectri_targets     = curated
)
print(tab[, c("regulator", "regulon_size_ARACNe", "z_stouffer", "transcript_replicated",
              "top_term", "collectri_targets")], row.names = FALSE)

write.table(tab, "results/integration/Table1_MR_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
