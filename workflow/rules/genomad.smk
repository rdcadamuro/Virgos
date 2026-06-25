# ======================================================================
# BLOCO 01 — geNomad: identifica contigs virais (genoma OU metagenoma)
# ======================================================================

rule genomad_end_to_end:
    input:
        fasta = s_assembly,
        db    = config["databases"]["genomad"],
    output:
        virus_fna = D_GENOMAD + "/{sample}/{sample}_summary/{sample}_virus.fna",
        virus_tsv = D_GENOMAD + "/{sample}/{sample}_summary/{sample}_virus_summary.tsv",
    params:
        outdir = lambda w: f"{D_GENOMAD}/{w.sample}",
        extra  = config["genomad"]["extra"],
    threads: config["threads"]["genomad"]
    conda: "../envs/genomad.yaml"
    log: D_LOGS + "/01_genomad/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {params.outdir}
        echo "=== geNomad | {wildcards.sample} | $(date) ===" > "$LOG"
        ln -sf $(readlink -f {input.fasta}) {params.outdir}/{wildcards.sample}.fna
        # geNomad sob set -e: se falhar (exit!=0), a regra falha e o Snakemake aborta.
        genomad end-to-end {params.extra} -t {threads} \
            {params.outdir}/{wildcards.sample}.fna {params.outdir} {input.db} \
            >> "$LOG" 2>&1
        # --- VERIFICACAO: chegou aqui => geNomad COMPLETOU (exit 0) ---
        nvir=$(grep -c '^>' {output.virus_fna} 2>/dev/null || echo 0)
        echo "[geNomad] OK -- rodou ate o fim (exit 0). Sequencias virais: $nvir" >> "$LOG"
        if [ "$nvir" -eq 0 ]; then
            echo "[geNomad] NOTA: 0 virus. geNomad COMPLETOU normalmente; a amostra nao tem sinal viral. NAO e erro." >> "$LOG"
        fi
        """


rule genomad_filter_score:
    input:
        fna = rules.genomad_end_to_end.output.virus_fna,
        tsv = rules.genomad_end_to_end.output.virus_tsv,
    output:
        fna = VIRAL_FNA,
        ids = D_GENOMAD + "/viral_contigs/{sample}.keep_ids.txt",
    params:
        min_score = config["genomad"]["min_virus_score"],
        script    = f"{config['scripts_dir']}/filter_genomad.py",
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        python {params.script} \
            --fna {input.fna} --summary {input.tsv} \
            --min-score {params.min_score} \
            --out-fna {output.fna} --out-ids {output.ids}
        """
