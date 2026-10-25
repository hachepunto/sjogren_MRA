# Reconstruye resources/collectri_regulons_human.tsv desde el servicio web de
# OmniPath (fuera del flujo de Snakemake: el archivo se distribuye con el repo).
# Complejos excluidos; mor = -1 si la interacción solo es inhibitoria, +1 en otro caso.
# Uso: Rscript scripts/get_collectri.R
url <- "https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=1"
x <- read.delim(url)
x <- x[!grepl("^COMPLEX:", x$source), ]
ct <- unique(data.frame(source = x$source_genesymbol, target = x$target_genesymbol,
                        mor = ifelse(x$is_inhibition == "True" & x$is_stimulation == "False", -1L, 1L)))
write.table(ct, "resources/collectri_regulons_human.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(ct), "interacciones,", length(unique(ct$source)), "TFs\n")
