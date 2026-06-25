# ======================================================================
# BLOCO 02 — CheckV: QC sobre os CONTIGS virais do geNomad.
# Roda 1x por amostra/genoma e produz viral_genes/completeness/quality
# por contig. O filtro de "tem hallmark (>=1 gene viral)" é aplicado no
# bloco 03 (em bins p/ metagenoma, em contigs p/ genoma).
# ======================================================================

rule checkv:
    input:
        fna = VIRAL_FNA,
        db  = config["databases"]["checkv"],
    output:
        summary = CHECKV_SUM,
    params:
        outdir = lambda w: f"{D_CHECKV}/{w.sample}",
    threads: config["threads"]["checkv"]
    conda: "../envs/checkv.yaml"
    log: D_LOGS + "/02_checkv/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {params.outdir}
        echo "=== CheckV | {wildcards.sample} | $(date) ===" > "$LOG"
        if [ ! -s {input.fna} ]; then
            echo "[CheckV] SEM contigs virais de entrada para {wildcards.sample} -> tabela vazia. NAO e erro (geNomad nao achou virus)." >> "$LOG"
            echo -e "contig_id\tcontig_length\tprovirus\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcontamination" > {output.summary}
            exit 0
        fi
        checkv end_to_end {input.fna} {params.outdir} \
            -t {threads} -d {input.db} >> "$LOG" 2>&1
        # --- VERIFICACAO: chegou aqui => CheckV COMPLETOU (exit 0) ---
        n=$(tail -n +2 {output.summary} 2>/dev/null | wc -l)
        echo "[CheckV] OK -- rodou ate o fim (exit 0). Contigs avaliados: $n" >> "$LOG"
        """


# ======================================================================
# CheckV no BIN INTEIRO (vMAG) — scaffolds do bin unidos por N's (link_bins).
# Recomendacao do artigo do vRhyme: estima completude/qualidade do bin como
# uma unidade, ignorando a fragmentacao. CheckV end_to_end chama prodigal-gv
# internamente; o espacador de N's impede genes atravessando as juncoes.
# Usado SO no caminho com reads (binning). Genoma sem reads continua por-contig.
# ======================================================================
rule checkv_bins:
    input:
        fna = LINKED_BINS,
        db  = config["databases"]["checkv"],
    output:
        summary = CHECKV_BINS_SUM,
    params:
        outdir = lambda w: f"{D_CHECKV}/{w.sample}_bins",
    threads: config["threads"]["checkv"]
    conda: "../envs/checkv.yaml"
    log: D_LOGS + "/02_checkv_bins/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {params.outdir}
        echo "=== CheckV (BINS unidos por N) | {wildcards.sample} | $(date) ===" > "$LOG"
        if [ ! -s {input.fna} ]; then
            echo "[CheckV-bins] SEM bins de entrada para {wildcards.sample} -> tabela vazia. NAO e erro (0 bins formados)." >> "$LOG"
            echo -e "contig_id\tcontig_length\tprovirus\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcontamination" > {output.summary}
            exit 0
        fi
        checkv end_to_end {input.fna} {params.outdir} \
            -t {threads} -d {input.db} >> "$LOG" 2>&1
        n=$(tail -n +2 {output.summary} 2>/dev/null | wc -l)
        echo "[CheckV-bins] OK -- rodou ate o fim (exit 0). Bins (vMAGs) avaliados: $n" >> "$LOG"
        """
