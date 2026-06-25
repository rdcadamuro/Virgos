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
