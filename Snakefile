configfile: "config.yaml"

R = config["rscript"]
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ARACNE_DIR = "results/networks/aracne"

wildcard_constraints:
    gse = r"GSE\d+",
    seed = r"\d+",


rule all:
    input:
        "results/figures/Top10_MRs_GSE51092.pdf",
        "results/figures/MRs_scatter_NES.pdf",
        "results/figures/MRs_integrative_figure.pdf",
        "results/figures/ORA_combined.pdf",
        "results/figures/DE_integrative_figure.pdf",
        "results/figures/ssGSEA_GSE140161_replicated.pdf",
        "results/integration/MR_DE_convergence.txt",
        "results/integration/meta_DE_results.rda",
        "results/validation/Table_S1_DE_GSE66795_replication.xlsx",
        "results/robustness/viper_cross_network.rda",
        "results/robustness/viper_collectri.rda",
        "results/robustness/deconvolution_mra.rda",
        "results/enrichment/fgsea_regulons.rda",


# ── Descargas ────────────────────────────────────────────────────────────────
rule series_matrix:
    output: "data/{gse}_series_matrix.txt"
    params: url = lambda w: f"{GEO}/{w.gse[:-3]}nnn/{w.gse}/matrix/{w.gse}_series_matrix.txt.gz"
    shell: "curl -sSfL {params.url} | gunzip -c > {output}"

rule gse140161_processed:
    output: "data/GSE140161_Processed_data.txt.gz"
    shell: "curl -sSfL -o {output} {GEO}/GSE140nnn/GSE140161/suppl/GSE140161_Processed_data.txt.gz"

rule lambert_csv:
    output: "data/lambert2018_DatabaseExtract_v1.01.csv"
    shell: "curl -sSfL -o {output} https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"


# ── Redes ────────────────────────────────────────────────────────────────────
rule tf_list:
    input: "data/lambert2018_DatabaseExtract_v1.01.csv"
    output: "results/networks/TFs_Lambert2018.txt"
    log: "logs/tf_list.log"
    shell: "{R} scripts/tf_list_lambert.R > {log} 2>&1"

rule prep_gse51092:
    input: "data/GSE51092_series_matrix.txt", "resources/probe_map_GSE51092_2019.tsv"
    output:
        eset = "results/mra/eset_GSE51092.rda",
        cases = "results/networks/SJS_matrix.txt",
    log: "logs/prep_gse51092.log"
    shell: "{R} scripts/prep_GSE51092.R > {log} 2>&1"

if config["run_aracne"]:
    ARACNE = f"java -Xmx{config['aracne_mem']} -jar {config['aracne_jar']}"

    rule aracne_threshold:
        input:
            exp = "results/networks/SJS_matrix.txt",
            tfs = "results/networks/TFs_Lambert2018.txt",
        output: directory(f"{ARACNE_DIR}/threshold")
        log: "logs/aracne/threshold.log"
        threads: 3
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
        resources: mem_mb = 4000
        shell:
            "rm -rf {params.tmp} && mkdir -p {params.tmp} && cp {input.thr}/miThreshold_* {params.tmp}/ && "
            "{ARACNE} -e {input.exp} -o {params.tmp} --tfs {input.tfs} "
            "--pvalue {config[aracne_pvalue]} --seed {wildcards.seed} > {log} 2>&1 && "
            "mv {params.tmp}/bootstrapNetwork_*.txt {output} && rm -rf {params.tmp}"

    rule aracne_consolidate:
        input: expand(f"{ARACNE_DIR}/bootstraps/bootstrapNetwork_seed{{seed}}.txt", seed=range(1, config["aracne_bootstraps"] + 1))
        output: "results/networks/sjs_network.txt"
        params: dir = f"{ARACNE_DIR}/bootstraps"
        log: "logs/aracne/consolidate.log"
        threads: 3
        shell:
            "{ARACNE} -o {params.dir} --consolidate > {log} 2>&1 && "
            "cp {params.dir}/network.txt {output}"
else:
    rule aracne_frozen:
        input: "resources/sjs_network_aracne.txt.gz"
        output: "results/networks/sjs_network.txt"
        shell: "gunzip -c {input} > {output}"


# ── MRA y expresión diferencial por cohorte ──────────────────────────────────
rule mra_gse51092:
    input:
        "results/mra/eset_GSE51092.rda",
        "results/networks/sjs_network.txt",
    output:
        "results/mra/mra_GSE51092.rda",
        "results/mra/top100_MRs_GSE51092.txt",
        "results/de/DE_limma_GSE51092.txt",
        "results/figures/Top10_MRs_GSE51092.pdf",
    log: "logs/mra_gse51092.log"
    shell: "{R} scripts/mra_GSE51092.R > {log} 2>&1"

rule mra_gse84844:
    input:
        "data/GSE84844_series_matrix.txt",
        "results/networks/TFs_Lambert2018.txt",
    output:
        "results/networks/regulon_GSE84844_corto.rda",
        "results/mra/mra_GSE84844.rda",
        "results/mra/top100_MRs_GSE84844.txt",
        "results/de/DE_limma_GSE84844.txt",
    threads: config["corto_threads"]
    log: "logs/mra_gse84844.log"
    # los hilos son los workers de corto; BLAS en 1 para que no se multipliquen
    shell: "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {R} scripts/mra_GSE84844.R {threads} > {log} 2>&1"


# ── Integración entre cohortes ───────────────────────────────────────────────
rule compare_mrs:
    input: "results/mra/mra_GSE51092.rda", "results/mra/mra_GSE84844.rda"
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
    output:
        "results/integration/MRs_meta_shared.txt",
        "results/integration/MRs_concordance_tests.txt",
        "results/integration/meta_MRs_results.rda",
    log: "logs/meta_mrs.log"
    shell: "{R} scripts/meta_MRs.R > {log} 2>&1"

rule figure_mrs:
    input: "results/integration/comparison_results.rda", "results/integration/MRs_meta_shared.txt"
    output: "results/figures/MRs_integrative_figure.pdf"
    log: "logs/figure_mrs.log"
    shell: "{R} scripts/figure_MRs.R > {log} 2>&1"

rule ora_regulons:
    input: "results/mra/mra_GSE51092.rda", "results/integration/MRs_meta_shared.txt"
    output:
        "results/enrichment/ORA_GO_BP.txt",
        "results/enrichment/ORA_Reactome.txt",
        "results/enrichment/ORA_results.rda",
        "results/figures/ORA_combined.pdf",
    log: "logs/ora_regulons.log"
    shell: "{R} scripts/ora_regulons.R > {log} 2>&1"

rule compare_de:
    input: "results/de/DE_limma_GSE51092.txt", "results/de/DE_limma_GSE84844.txt"
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
    input: "results/de/comparison_DE_results.rda"
    output: "results/figures/DE_integrative_figure.pdf"
    log: "logs/figure_de.log"
    shell: "{R} scripts/figure_DE.R > {log} 2>&1"

rule meta_de:
    input:
        "results/de/DE_comparacion.txt",
        "results/de/DE_RRA_upregulated.txt",
        "results/de/DE_RRA_downregulated.txt",
    output:
        "results/integration/DE_meta_shared.txt",
        "results/integration/meta_DE_results.rda",
    log: "logs/meta_de.log"
    shell: "{R} scripts/meta_DE.R > {log} 2>&1"

rule convergence:
    input: "results/integration/MRs_meta_shared.txt", "results/integration/DE_meta_shared.txt"
    output:
        "results/integration/MR_DE_convergence.txt",
        "results/integration/MR_DE_convergence_summary.txt",
    log: "logs/convergence.log"
    shell: "{R} scripts/convergence_MR_DE.R > {log} 2>&1"


# ── Cohortes de validación ───────────────────────────────────────────────────
rule de_gse66795:
    input: "data/GSE66795_series_matrix.txt", "results/de/comparison_DE_results.rda"
    output:
        "results/validation/DE_limma_GSE66795.txt",
        "results/validation/DE_GSE66795_replication.rda",
        "results/validation/Table_S1_DE_GSE66795_replication.xlsx",
    log: "logs/de_gse66795.log"
    shell: "{R} scripts/de_GSE66795.R > {log} 2>&1"

rule ssgsea_gse140161:
    input: "data/GSE140161_Processed_data.txt.gz", "results/de/DE_replicated_up.txt"
    output:
        "results/validation/ssGSEA_GSE140161_scores_replicated.txt",
        "results/validation/ssGSEA_GSE140161_replicated.rda",
        "results/figures/ssGSEA_GSE140161_replicated.pdf",
    log: "logs/ssgsea_gse140161.log"
    shell: "{R} scripts/ssgsea_GSE140161.R > {log} 2>&1"


# ── Análisis de robustez ─────────────────────────────────────────────────────
MRA = ["results/mra/mra_GSE51092.rda", "results/mra/mra_GSE84844.rda"]
CORTO = "results/networks/regulon_GSE84844_corto.rda"
CONSENSUS = "results/integration/MRs_meta_shared.txt"

rule viper_cross_network:
    input: MRA, CORTO, CONSENSUS
    output:
        "results/robustness/viper_cross_network.txt",
        "results/robustness/viper_cross_network_concordance.txt",
        "results/robustness/viper_cross_network.rda",
    log: "logs/viper_cross_network.log"
    shell: "{R} scripts/viper_cross_network.R > {log} 2>&1"

rule viper_collectri:
    input: MRA, CONSENSUS, "resources/collectri_regulons_human.tsv"
    output:
        "results/robustness/viper_collectri.txt",
        "results/robustness/viper_collectri_concordance.txt",
        "results/robustness/consensus_vs_collectri.txt",
        "results/robustness/viper_collectri.rda",
    log: "logs/viper_collectri.log"
    shell: "{R} scripts/viper_collectri.R > {log} 2>&1"

rule deconvolution:
    input: MRA, CORTO, CONSENSUS
    output:
        "results/robustness/deconv_composition_tests.txt",
        "results/robustness/mra_adjusted_composition.txt",
        "results/robustness/mra_adjusted_concordance.txt",
        "results/robustness/deconvolution_mra.rda",
    log: "logs/deconvolution.log"
    shell: "{R} scripts/deconvolution_mra.R > {log} 2>&1"

rule fgsea_regulons:
    input: MRA, CORTO, CONSENSUS
    output:
        "results/enrichment/fgsea_regulons_long.txt",
        "results/enrichment/fgsea_regulons_wide.txt",
        "results/enrichment/fgsea_vs_viper_concordance.txt",
        "results/enrichment/fgsea_regulons.rda",
    log: "logs/fgsea_regulons.log"
    shell: "{R} scripts/fgsea_regulons.R > {log} 2>&1"
