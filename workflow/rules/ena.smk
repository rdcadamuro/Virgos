# ======================================================================
# BLOCO ENA — pacote de submissao dos uViGs ao European Nucleotide Archive.
#
# Auto-alimenta as tabelas MIUViG (Minimum Information about an Uncultivated
# Virus Genome) + manifestos Webin-CLI (assembly tipo MAG) + FASTAs + guia,
# a partir do que o pipeline JA produz (CheckV, VITAP, derep, geNomad, mobiloma).
# Modular e mode-aware (genoma ou metagenoma). Degrada com elegancia.
#
# Saidas em {D_SUM}/ena_submission/:
#   ena_uvig_metadata.tsv  ena_sample_registration.tsv  ena_samples.xml
#   manifests/*.manifest    fasta/*.fasta.gz   submit_all.sh
#   ENA_SUBMISSION_GUIDE.md
# ======================================================================

_ENA = config.get("ena", {})
D_ENA = D_SUM + "/ena_submission"


# Serializa a config 'ena' (+ parametros derivados de outras secoes) num JSON.
# Feito num 'run:' (Python no processo do Snakemake) p/ NAO injetar chaves '{}'
# de JSON na string de shell (quebraria o str.format do Snakemake).
rule ena_config:
    output:
        cfg = D_ENA + "/ena_config.json",
    run:
        import json
        import os
        os.makedirs(os.path.dirname(output.cfg), exist_ok=True)
        ena = dict(_ENA)
        ena.setdefault("genomad_min_virus_score",
                       config.get("genomad", {}).get("min_virus_score", 0.7))
        ena.setdefault("derep_ani", config.get("derep", {}).get("min_ani", 95))
        ena.setdefault("derep_af", config.get("derep", {}).get("min_af", 85))
        with open(output.cfg, "w") as fh:
            json.dump(ena, fh, indent=2)


def _ena_abund_arg(wildcards):
    """Em metagenoma, passa a abundancia (CoverM) p/ COVERAGE; senao vazio
    (o script cai p/ o cov no nome do contig)."""
    if MODE == "metagenome":
        return f"--abundance {D_ABUND}/votu_abundance_tpm.tsv"
    return ""


rule ena_tables:
    input:
        cfg         = rules.ena_config.output.cfg,
        per_bin     = D_SUM + "/viral_bins_table.tsv",
        checkv_bins = expand(CHECKV_BINS_SUM, sample=SAMPLES),
        vitap       = D_VITAP + "/vitap_taxonomy.tsv",
        votus       = VOTUS,
    output:
        meta  = D_ENA + "/ena_uvig_metadata.tsv",
        guide = D_ENA + "/ENA_SUBMISSION_GUIDE.md",
        xml   = D_ENA + "/ena_samples.xml",
    params:
        script   = f"{config['scripts_dir']}/build_ena_tables.py",
        outdir   = D_ENA,
        mode     = MODE,
        mobilome = D_ANNOT + "/mobilome_table.tsv",
        samples  = config["samples"],
        abund    = _ena_abund_arg,
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/09_ena.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {params.outdir}
        echo "=== ENA submission package | $(date) ===" > "$LOG"
        # mobiloma e opcional (is_provirus tambem vem do CheckV/nome) -> so passa se existir
        MOBARG=""
        [ -s "{params.mobilome}" ] && MOBARG="--mobilome {params.mobilome}"
        python {params.script} \
            --per-bin {input.per_bin} \
            --checkv-bins {input.checkv_bins} \
            --vitap {input.vitap} \
            --votus {input.votus} \
            --samples {params.samples} \
            --config-json {input.cfg} \
            --mode {params.mode} \
            $MOBARG {params.abund} \
            --outdir {params.outdir} >> "$LOG" 2>&1
        n=$(tail -n +2 {output.meta} 2>/dev/null | wc -l)
        m=$(ls {params.outdir}/manifests/*.manifest 2>/dev/null | wc -l)
        echo "[ENA] OK -- rodou ate o fim (exit 0). uViGs: $n; manifestos: $m" >> "$LOG"
        """
