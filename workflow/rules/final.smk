# ======================================================================
# BLOCO FINAL — Tabela por vOTU + resumo por amostra.
# metagenome -> 07_summary_final (inclui prevalência/abundância dos vOTUs)
# genome     -> 06_summary_final
# ======================================================================

_MODE = config.get("mode", "metagenome")


def final_inputs(wildcards):
    d = dict(
        clusters = CLUSTERS,
        votu_map = VOTU_MAP,
        units    = ALL_UNITS,
        phatyp   = D_PHABOX + "/final_prediction/phatyp_prediction.tsv",
        vitap    = D_VITAP + "/vitap_taxonomy.tsv",
    )
    if _MODE == "metagenome":
        d["abund_tpm"] = D_ABUND + "/votu_abundance_tpm.tsv"
    return d


def abund_arg(wildcards, input):
    tpm = getattr(input, "abund_tpm", None)
    return f"--abundance-tpm {tpm}" if tpm else ""


rule final_table:
    input:
        unpack(final_inputs)
    output:
        per_bin    = D_SUM + "/viral_bins_table.tsv",
        per_sample = D_SUM + "/per_sample_summary.tsv",
    params:
        script    = f"{config['scripts_dir']}/build_final_table.py",
        abund_arg = abund_arg,
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        python {params.script} \
            --clusters {input.clusters} --votu-map {input.votu_map} \
            --units {input.units} --phatyp {input.phatyp} --vitap {input.vitap} \
            {params.abund_arg} \
            --out-per-bin {output.per_bin} --out-per-sample {output.per_sample}
        """


# ----------------------------------------------------------------------
# CANDIDATO A ESPÉCIE NOVA — flag honesto/conservador por vOTU.
# Só marca 'candidate_novel_species' se High-quality+ (a def. de espécie
# 95%ANI/85%AF só é confiável fora de fragmento) + fago + sem táxon de
# gênero/espécie no VITAP. Fragmento -> 'indeterminate_fragment'.
# ----------------------------------------------------------------------
rule votu_novelty:
    input:
        per_bin     = D_SUM + "/viral_bins_table.tsv",
        checkv_bins = expand(CHECKV_BINS_SUM, sample=SAMPLES),
    output:
        tsv = D_SUM + "/votu_novelty.tsv",
    params:
        script = f"{config['scripts_dir']}/flag_novelty.py",
        minc   = config.get("annotation", {}).get("novelty_min_completeness", 90),
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/08_novelty.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.tsv})"
        echo "=== candidato a especie nova | $(date) ===" > "$LOG"
        python {params.script} --per-bin {input.per_bin} \
            --checkv-bins {input.checkv_bins} --min-completeness {params.minc} \
            --out {output.tsv} >> "$LOG" 2>&1
        n=$(grep -c 'candidate_novel_species' {output.tsv} 2>/dev/null || echo 0)
        echo "[novelty] OK -- rodou ate o fim (exit 0). candidatos a especie nova: $n" >> "$LOG"
        """


# ======================================================================
# VERIFICACAO INTELIGENTE — roda DEPOIS de tudo. Audita todas as etapas,
# cruza consistencia entre elas e varre os logs por erros silenciosos.
# Sempre completa (nao-strict): grava RUN_HEALTHCHECK.txt + imprime no console.
# ======================================================================
rule healthcheck:
    input:
        per_bin    = D_SUM + "/viral_bins_table.tsv",
        per_sample = D_SUM + "/per_sample_summary.tsv",
    output:
        report = D_SUM + "/RUN_HEALTHCHECK.txt",
    params:
        script  = f"{config['scripts_dir']}/verify_run.py",
        outdir  = OUT,
        mode    = _MODE,
        samples = config["samples"],
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        python {params.script} --outdir {params.outdir} --mode {params.mode} \
            --samples {params.samples} --out {output.report}
        """
