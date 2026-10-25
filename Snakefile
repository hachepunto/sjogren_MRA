configfile: "config.yaml"

R = config["rscript"]
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
ARACNE_DIR = "results/networks/aracne"

# BLAS, OpenMP, readr y data.table en un hilo: cada regla declara en threads:
# las CPUs que de verdad usa
shell.prefix("set -euo pipefail; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "
             "VROOM_THREADS=1 R_DATATABLE_NUM_THREADS=1; ")

wildcard_constraints:
    gse = r"GSE\d+",
    gpl = r"GPL\d+",
    seed = r"\d+",


def code(name, helpers=False):
    """Script de la regla (y helpers.R) como entrada: un cambio de código la rehace."""
    return [f"scripts/{name}.R"] + (["scripts/helpers.R"] if helpers else [])


rule all:
    input:
        "results/data_checksums.txt",
        "results/figures/Top10_MRs_GSE51092.pdf",
        "results/figures/MRs_scatter_NES.pdf",
        "results/figures/MRs_integrative_figure.pdf",
        "results/figures/ORA_combined.pdf",
        "results/figures/DE_integrative_figure.pdf",
        "results/figures/ssGSEA_GSE140161_replicated.pdf",
        "results/integration/MR_DE_convergence.txt",
        "results/integration/meta_DE_results.rda",
        "results/integration/Table1_MR_summary.tsv",
        "results/validation/Table_S1_DE_GSE66795_replication.xlsx",
        "results/enrichment/ORA_nonconsensus.tsv",
        "results/robustness/composition_reported.tsv",
        "results/Supplementary_Tables_S1-S7.xlsx",


# ── Descargas ────────────────────────────────────────────────────────────────
rule series_matrix:
    output: "data/{gse}_series_matrix.txt"
    params: url = lambda w: f"{GEO}/series/{w.gse[:-3]}nnn/{w.gse}/matrix/{w.gse}_series_matrix.txt.gz"
    log: "logs/download/{gse}_series_matrix.log"
    shell: "curl -sSfL {params.url} 2> {log} | gunzip -c > {output}"

# Anotación de plataforma con el nombre y formato que GEOquery::getGEO() busca en
# data/: la versión completa (texto plano, aunque termine en .gz) o la anotación
# curada .annot.gz
rule gpl_soft:
    output: "data/{gpl}.soft.gz"
    params: url = lambda w: f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={w.gpl}&form=text&view=full"
    log: "logs/download/{gpl}_soft.log"
    shell: "curl -sSfL -o {output} '{params.url}' 2> {log}"

rule gpl_annot:
    output: "data/{gpl}.annot.gz"
    params: url = lambda w: f"{GEO}/platforms/{w.gpl[:-3]}nnn/{w.gpl}/annot/{w.gpl}.annot.gz"
    log: "logs/download/{gpl}_annot.log"
    shell: "curl -sSfL -o {output} {params.url} 2> {log}"

rule gse140161_processed:
    output: "data/GSE140161_Processed_data.txt.gz"
    log: "logs/download/GSE140161_processed.log"
    shell: "curl -sSfL -o {output} {GEO}/series/GSE140nnn/GSE140161/suppl/GSE140161_Processed_data.txt.gz 2> {log}"

rule lambert_csv:
    output: "data/lambert2018_DatabaseExtract_v1.01.csv"
    log: "logs/download/lambert_csv.log"
    shell: "curl -sSfL -o {output} https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv 2> {log}"

# Firma de MCP-counter del commit instalado (envs/install_r_packages.R)
rule mcpcounter_genes:
    output: "data/MCPcounter_genes.txt"
    log: "logs/download/mcpcounter_genes.log"
    shell:
        "curl -sSfL -o {output} https://raw.githubusercontent.com/ebecht/MCPcounter/"
        "b6eac73e91c246fcff0bb1a5c68a816cd588fc48/Signatures/genes.txt 2> {log}"

DATA = [
    "data/GSE51092_series_matrix.txt", "data/GSE84844_series_matrix.txt",
    "data/GSE66795_series_matrix.txt", "data/GSE140161_Processed_data.txt.gz",
    "data/GPL6884.soft.gz", "data/GPL570.annot.gz", "data/GPL10558.soft.gz",
    "data/lambert2018_DatabaseExtract_v1.01.csv", "data/MCPcounter_genes.txt",
]

# Compara las descargas con las del manuscrito; si GEO actualizó un archivo lo
# informa en results/data_checksums.txt sin detener el flujo
rule data_checksums:
    input: DATA, md5 = "resources/data_md5.txt"
    output: "results/data_checksums.txt"
    log: "logs/data_checksums.log"
    shell:
        "md5sum -c {input.md5} > {output} 2>&1 || "
        "echo 'WARNING: some downloads differ from the files used for the manuscript' | tee -a {output} > {log}"


# ── Redes ────────────────────────────────────────────────────────────────────
rule tf_list:
    input: "data/lambert2018_DatabaseExtract_v1.01.csv", code("tf_list_lambert")
    output: "results/networks/TFs_Lambert2018.txt"
    log: "logs/tf_list.log"
    shell: "{R} scripts/tf_list_lambert.R > {log} 2>&1"

rule prep_gse51092:
    input:
        "data/GSE51092_series_matrix.txt", "data/GPL6884.soft.gz",
        "resources/probe_map_GSE51092_2019.tsv", code("prep_GSE51092", helpers=True),
    output:
        eset = "results/mra/eset_GSE51092.rda",
        cases = "results/networks/SJS_matrix.txt",
    log: "logs/prep_gse51092.log"
    shell: "{R} scripts/prep_GSE51092.R > {log} 2>&1"

if config["run_aracne"]:
    ARACNE = f"java -Xmx{config['aracne_mem']} -jar {config['aracne_jar']}"
    # RAM por JVM: el heap de aracne_mem más ~25 % de memoria fuera del heap
    ARACNE_MEM_MB = int(float(config["aracne_mem"].rstrip("GgMm")) *
                        (1024 if config["aracne_mem"][-1] in "Gg" else 1) * 1.25)
    N_BOOT = config["aracne_bootstraps"]

    rule aracne_threshold:
        input:
            exp = "results/networks/SJS_matrix.txt",
            tfs = "results/networks/TFs_Lambert2018.txt",
        output: directory(f"{ARACNE_DIR}/threshold")
        log: "logs/aracne/threshold.log"
        threads: 3
        resources: mem_mb = ARACNE_MEM_MB
        shell:
            "{ARACNE} -e {input.exp} -o {output} --tfs {input.tfs} "
            "--pvalue {config[aracne_pvalue]} --seed 1 --calculateThreshold > {log} 2>&1"

    # Cada bootstrap corre en su propia carpeta (con copia del umbral) y su
    # salida, de nombre aleatorio, se renombra por semilla
    rule aracne_bootstrap:
        input:
            exp = "results/networks/SJS_matrix.txt",
            tfs = "results/networks/TFs_Lambert2018.txt",
            thr = f"{ARACNE_DIR}/threshold",
        output: f"{ARACNE_DIR}/bootstraps/bootstrapNetwork_seed{{seed}}.txt"
        params: tmp = f"{ARACNE_DIR}/tmp_seed{{seed}}"
        log: "logs/aracne/bootstrap_{seed}.log"
        threads: 3   # medido: ~2.2 CPUs por JVM (hilos de GC y JIT)
        resources: mem_mb = ARACNE_MEM_MB
        shell:
            "rm -rf {params.tmp} && mkdir -p {params.tmp} && cp {input.thr}/miThreshold_* {params.tmp}/ && "
            "{ARACNE} -e {input.exp} -o {params.tmp} --tfs {input.tfs} "
            "--pvalue {config[aracne_pvalue]} --seed {wildcards.seed} > {log} 2>&1 && "
            "mv {params.tmp}/bootstrapNetwork_*.txt {output} && rm -rf {params.tmp}"

    # ARACNe-AP consolida todos los bootstraps de la carpeta: se comprueba que
    # sean exactamente los N_BOOT esperados (sin restos de corridas anteriores)
    rule aracne_consolidate:
        input: expand(f"{ARACNE_DIR}/bootstraps/bootstrapNetwork_seed{{seed}}.txt", seed=range(1, N_BOOT + 1))
        output: "results/networks/sjs_network.txt"
        params: dir = f"{ARACNE_DIR}/bootstraps"
        log: "logs/aracne/consolidate.log"
        threads: 3
        resources: mem_mb = ARACNE_MEM_MB
        shell:
            "{{ test $(ls {params.dir}/bootstrapNetwork_*.txt | wc -l) -eq {N_BOOT} || "
            "{{ echo 'unexpected bootstrap files in {params.dir}' > {log}; exit 1; }}; }} && "
            "{ARACNE} -o {params.dir} --consolidate > {log} 2>&1 && "
            "mv {params.dir}/network.txt {output}"
else:
    rule aracne_frozen:
        input: "resources/sjs_network_aracne.txt.gz"
        output: "results/networks/sjs_network.txt"
        log: "logs/aracne_frozen.log"
        shell: "gunzip -c {input} > {output} 2> {log}"


# ── MRA y expresión diferencial por cohorte ──────────────────────────────────
rule mra_gse51092:
    input:
        "results/mra/eset_GSE51092.rda",
        "results/networks/sjs_network.txt",
        code("mra_GSE51092", helpers=True),
    output:
        "results/mra/mra_GSE51092.rda",
        "results/mra/top100_MRs_GSE51092.txt",
        "results/de/DE_limma_GSE51092.txt",
        "results/figures/Top10_MRs_GSE51092.pdf",
    log: "logs/mra_gse51092.log"
    shell: "{R} scripts/mra_GSE51092.R > {log} 2>&1"

rule mra_gse84844:
    input:
        "data/GSE84844_series_matrix.txt", "data/GPL570.annot.gz",
        "results/networks/TFs_Lambert2018.txt",
        code("mra_GSE84844", helpers=True),
    output:
        "results/networks/regulon_GSE84844_corto.rda",
        "results/mra/mra_GSE84844.rda",
        "results/mra/top100_MRs_GSE84844.txt",
        "results/de/DE_limma_GSE84844.txt",
    threads: config["corto_threads"]   # workers de corto
    log: "logs/mra_gse84844.log"
    shell: "{R} scripts/mra_GSE84844.R {threads} > {log} 2>&1"


# ── Integración entre cohortes ───────────────────────────────────────────────
MRA = ["results/mra/mra_GSE51092.rda", "results/mra/mra_GSE84844.rda"]
CORTO = "results/networks/regulon_GSE84844_corto.rda"
MRS_META = "results/integration/MRs_meta_shared.txt"
DE_META = "results/integration/DE_meta_shared.txt"

rule compare_mrs:
    input: MRA, code("compare_MRs")
    output:
        "results/integration/MRs_comparacion.txt",
        "results/integration/MRs_RRA_activated.txt",
        "results/integration/MRs_RRA_repressed.txt",
        "results/integration/comparison_results.rda",
        "results/figures/MRs_scatter_NES.pdf",
    log: "logs/compare_mrs.log"
    shell: "{R} scripts/compare_MRs.R > {log} 2>&1"

rule meta_mrs:
    input:
        "results/integration/MRs_comparacion.txt",
        "results/integration/MRs_RRA_activated.txt",
        "results/integration/MRs_RRA_repressed.txt",
        code("meta_MRs"),
    output:
        MRS_META,
        "results/integration/MRs_concordance_tests.txt",
        "results/integration/meta_MRs_results.rda",
    log: "logs/meta_mrs.log"
    shell: "{R} scripts/meta_MRs.R > {log} 2>&1"

rule figure_mrs:
    input: "results/integration/comparison_results.rda", MRS_META, code("figure_MRs")
    output: "results/figures/MRs_integrative_figure.pdf"
    log: "logs/figure_mrs.log"
    shell: "{R} scripts/figure_MRs.R > {log} 2>&1"

rule ora_regulons:
    input: "results/mra/mra_GSE51092.rda", MRS_META, code("ora_regulons", helpers=True)
    output:
        "results/enrichment/ORA_GO_BP.txt",
        "results/enrichment/ORA_Reactome.txt",
        "results/enrichment/ORA_results.rda",
        "results/figures/ORA_combined.pdf",
    log: "logs/ora_regulons.log"
    shell: "{R} scripts/ora_regulons.R > {log} 2>&1"

rule ora_nonconsensus:
    input: "results/mra/mra_GSE51092.rda", code("ora_nonconsensus", helpers=True)
    output:
        "results/enrichment/ORA_nonconsensus.tsv",
        "results/enrichment/ORA_nonconsensus.rda",
    log: "logs/ora_nonconsensus.log"
    shell: "{R} scripts/ora_nonconsensus.R > {log} 2>&1"

rule compare_de:
    input: "results/de/DE_limma_GSE51092.txt", "results/de/DE_limma_GSE84844.txt", code("compare_ED", helpers=True)
    output:
        "results/de/DE_comparacion.txt",
        "results/de/DE_RRA_upregulated.txt",
        "results/de/DE_RRA_downregulated.txt",
        "results/de/DE_replicated_up.txt",
        "results/de/DE_replicated_down.txt",
        "results/de/comparison_DE_results.rda",
        "results/enrichment/DE_ORA_GO_BP.txt",
        "results/enrichment/DE_ORA_Reactome.txt",
    log: "logs/compare_de.log"
    shell: "{R} scripts/compare_ED.R > {log} 2>&1"

rule figure_de:
    input: "results/de/comparison_DE_results.rda", code("figure_DE")
    output: "results/figures/DE_integrative_figure.pdf"
    log: "logs/figure_de.log"
    shell: "{R} scripts/figure_DE.R > {log} 2>&1"

rule meta_de:
    input:
        "results/de/DE_comparacion.txt",
        "results/de/DE_RRA_upregulated.txt",
        "results/de/DE_RRA_downregulated.txt",
        code("meta_DE"),
    output:
        DE_META,
        "results/integration/meta_DE_results.rda",
    log: "logs/meta_de.log"
    shell: "{R} scripts/meta_DE.R > {log} 2>&1"

rule convergence:
    input: MRS_META, DE_META, code("convergence_MR_DE")
    output:
        "results/integration/MR_DE_convergence.txt",
        "results/integration/MR_DE_convergence_summary.txt",
    log: "logs/convergence.log"
    shell: "{R} scripts/convergence_MR_DE.R > {log} 2>&1"


# ── Cohortes de validación ───────────────────────────────────────────────────
rule de_gse66795:
    input:
        "data/GSE66795_series_matrix.txt", "data/GPL10558.soft.gz",
        "results/de/comparison_DE_results.rda",
        code("de_GSE66795", helpers=True),
    output:
        "results/validation/DE_limma_GSE66795.txt",
        "results/validation/DE_GSE66795_replication.rda",
        "results/validation/Table_S1_DE_GSE66795_replication.xlsx",
    log: "logs/de_gse66795.log"
    shell: "{R} scripts/de_GSE66795.R > {log} 2>&1"

rule ssgsea_gse140161:
    input: "data/GSE140161_Processed_data.txt.gz", "results/de/DE_replicated_up.txt", code("ssgsea_GSE140161")
    output:
        "results/validation/ssGSEA_GSE140161_scores_replicated.txt",
        "results/validation/ssGSEA_GSE140161_replicated.rda",
        "results/figures/ssGSEA_GSE140161_replicated.pdf",
    log: "logs/ssgsea_gse140161.log"
    shell: "{R} scripts/ssgsea_GSE140161.R > {log} 2>&1"


# ── Análisis de robustez ─────────────────────────────────────────────────────
rule viper_cross_network:
    input: MRA, CORTO, MRS_META, code("viper_cross_network", helpers=True)
    output:
        "results/robustness/viper_cross_network.txt",
        "results/robustness/viper_cross_network_concordance.txt",
        "results/robustness/viper_cross_network.rda",
    log: "logs/viper_cross_network.log"
    shell: "{R} scripts/viper_cross_network.R > {log} 2>&1"

rule viper_collectri:
    input: MRA, MRS_META, "resources/collectri_regulons_human.tsv", code("viper_collectri", helpers=True)
    output:
        "results/robustness/viper_collectri.txt",
        "results/robustness/viper_collectri_concordance.txt",
        "results/robustness/consensus_vs_collectri.txt",
        "results/robustness/viper_collectri.rda",
    log: "logs/viper_collectri.log"
    shell: "{R} scripts/viper_collectri.R > {log} 2>&1"

rule deconvolution:
    input: MRA, CORTO, MRS_META, "data/MCPcounter_genes.txt", code("deconvolution_mra", helpers=True)
    output:
        "results/robustness/deconv_composition_tests.txt",
        "results/robustness/mra_adjusted_composition.txt",
        "results/robustness/mra_adjusted_concordance.txt",
        expand("results/robustness/{method}_{gse}.txt", method=["xcell", "mcpcounter"], gse=["GSE51092", "GSE84844"]),
        "results/robustness/deconvolution_mra.rda",
    log: "logs/deconvolution.log"
    shell: "{R} scripts/deconvolution_mra.R > {log} 2>&1"

rule composition_reported:
    input: "results/robustness/deconv_composition_tests.txt", code("composition_reported")
    output: "results/robustness/composition_reported.tsv"
    log: "logs/composition_reported.log"
    shell: "{R} scripts/composition_reported.R > {log} 2>&1"

rule fgsea_regulons:
    input: MRA, CORTO, MRS_META, code("fgsea_regulons", helpers=True)
    output:
        "results/enrichment/fgsea_regulons_long.txt",
        "results/enrichment/fgsea_regulons_wide.txt",
        "results/enrichment/fgsea_vs_viper_concordance.txt",
        "results/enrichment/fgsea_regulons.rda",
    log: "logs/fgsea_regulons.log"
    shell: "{R} scripts/fgsea_regulons.R > {log} 2>&1"


# ── Tablas del artículo ──────────────────────────────────────────────────────
rule table1_mr_summary:
    input:
        "results/mra/mra_GSE51092.rda", MRS_META, DE_META,
        "results/enrichment/ORA_GO_BP.txt",
        "results/enrichment/ORA_Reactome.txt",
        "results/robustness/consensus_vs_collectri.txt",
        code("table1_mr_summary"),
    output: "results/integration/Table1_MR_summary.tsv"
    log: "logs/table1_mr_summary.log"
    shell: "{R} scripts/table1_mr_summary.R > {log} 2>&1"

rule supplementary_tables:
    input:
        "results/validation/Table_S1_DE_GSE66795_replication.xlsx",
        MRS_META, DE_META,
        "results/robustness/viper_cross_network.txt",
        "results/robustness/viper_cross_network_concordance.txt",
        "results/robustness/viper_collectri.txt",
        "results/robustness/consensus_vs_collectri.txt",
        "results/robustness/viper_collectri_concordance.txt",
        expand("results/robustness/{method}_{gse}.txt", method=["xcell", "mcpcounter"], gse=["GSE51092", "GSE84844"]),
        "results/robustness/deconv_composition_tests.txt",
        "results/robustness/mra_adjusted_composition.txt",
        "results/robustness/mra_adjusted_concordance.txt",
        "results/enrichment/fgsea_regulons_long.txt",
        "results/enrichment/fgsea_vs_viper_concordance.txt",
        code("supplementary_tables"),
    output: "results/Supplementary_Tables_S1-S7.xlsx"
    log: "logs/supplementary_tables.log"
    shell: "{R} scripts/supplementary_tables.R > {log} 2>&1"
