# ======================================================================
# BLOCO 06 — Abundância de vOTUs (SÓ metagenoma).
# Mapeia os reads de cada amostra de volta nos vOTUs finais (já filtrados
# por CheckV) e calcula TPM, profundidade média e fração coberta.
# -> matriz vOTU x amostra (ecologia).
# ======================================================================

VOTU_BAM = D_ABUND + "/bam/{sample}.votu.sorted.bam"


rule map_to_votus:
    """Mapeia reads -> vOTUs (mesmo esquema de read_type do binning)."""
    input:
        ref = VOTUS,
        r1  = s_reads1,
    output:
        bam = VOTU_BAM,
        bai = VOTU_BAM + ".bai",
    params:
        read_type = s_readtype,
        r2        = s_reads2,
        mapper    = config["mapping"]["illumina_mapper"],
        preset    = config["mapping"]["nanopore_preset"],
        tmp       = lambda w: f"{D_ABUND}/bam/{w.sample}.tmp",
    threads: config["threads"]["mapping"]
    conda: "../envs/mapping.yaml"
    log: D_LOGS + "/06_map_votus/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam}) {params.tmp}
        RT="{params.read_type}"; REF="{input.ref}"
        if [ ! -s "$REF" ]; then : > {output.bam}; touch {output.bai}; exit 0; fi
        if [ "$RT" = "nano" ]; then
            minimap2 -t {threads} -ax {params.preset} "$REF" {input.r1} \
              | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
        elif [ "$RT" = "ilu-paired" ]; then
            if [ "{params.mapper}" = "bowtie2" ]; then
                bowtie2-build --threads {threads} "$REF" {params.tmp}/idx
                bowtie2 -p {threads} -x {params.tmp}/idx -1 {input.r1} -2 {params.r2} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            else
                bwa-mem2 index -p {params.tmp}/idx "$REF"
                bwa-mem2 mem -t {threads} {params.tmp}/idx {input.r1} {params.r2} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            fi
        else  # ilu-single
            if [ "{params.mapper}" = "bowtie2" ]; then
                bowtie2-build --threads {threads} "$REF" {params.tmp}/idx
                bowtie2 -p {threads} -x {params.tmp}/idx -U {input.r1} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            else
                bwa-mem2 index -p {params.tmp}/idx "$REF"
                bwa-mem2 mem -t {threads} {params.tmp}/idx {input.r1} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            fi
        fi
        samtools index -@ {threads} {output.bam}
        rm -rf {params.tmp}
        """


rule abundance_matrix:
    """CoverM agrega todos os BAMs num único passo -> matriz vOTU x amostra,
    para 3 métricas: TPM, profundidade média e fração coberta."""
    input:
        bams  = expand(VOTU_BAM, sample=SAMPLES),
        votus = VOTUS,
    output:
        tpm   = D_ABUND + "/votu_abundance_tpm.tsv",
        mean  = D_ABUND + "/votu_abundance_meandepth.tsv",
        covf  = D_ABUND + "/votu_abundance_coveredfraction.tsv",
    params:
        script = f"{config['scripts_dir']}/clean_coverm.py",
        # min_id cai p/ 85 automaticamente se houver amostra Nanopore (ver Snakefile)
        min_id = abundance_min_id(),
        min_aln = config.get("abundance", {}).get("min_read_aln_pct", 75),
    threads: config["threads"]["mapping"]
    conda: "../envs/abundance.yaml"
    log: D_LOGS + "/06_coverm.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tpm})
        if [ ! -s {input.votus} ]; then
            for f in {output.tpm} {output.mean} {output.covf}; do echo -e "votu" > $f; done
            exit 0
        fi
        for m in tpm mean covered_fraction; do
            coverm contig --bam-files {input.bams} \
                --min-read-percent-identity {params.min_id} \
                --min-read-aligned-percent {params.min_aln} \
                --methods $m --threads {threads} > {D_ABUND}/.coverm_$m.raw 2>> {log}
        done
        python {params.script} --raw {D_ABUND}/.coverm_tpm.raw --out {output.tpm}
        python {params.script} --raw {D_ABUND}/.coverm_mean.raw --out {output.mean}
        python {params.script} --raw {D_ABUND}/.coverm_covered_fraction.raw --out {output.covf}
        rm -f {D_ABUND}/.coverm_*.raw
        """
