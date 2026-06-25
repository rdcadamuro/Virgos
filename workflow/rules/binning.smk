# ======================================================================
# BLOCO 03 — Unidades virais finais (passed.fna).
# Binning é dirigido pela PRESENÇA DE READS (vale nos dois modos):
#   tem reads  -> mapeia -> vRhyme + VAMB -> bins -> filtra bins SEM hallmark
#   sem reads  -> usa os contigs virais c/ >=1 gene viral como unidades
# ======================================================================

# --------------------------- mapeamento -------------------------------
rule map_reads:
    input:
        ref = VIRAL_FNA,
        r1  = s_reads1,
    output:
        bam = BAM,
        bai = BAM + ".bai",
    params:
        read_type = s_readtype,
        r2        = s_reads2,
        mapper    = config["mapping"]["illumina_mapper"],
        preset    = config["mapping"]["nanopore_preset"],
        tmp       = lambda w: f"{D_BINS}/mapping/{w.sample}.tmp",
    threads: config["threads"]["mapping"]
    conda: "../envs/mapping.yaml"
    log: D_LOGS + "/03_mapping/{sample}.log"
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
        elif [ "$RT" = "ilu-single" ]; then
            if [ "{params.mapper}" = "bowtie2" ]; then
                bowtie2-build --threads {threads} "$REF" {params.tmp}/idx
                bowtie2 -p {threads} -x {params.tmp}/idx -U {input.r1} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            else
                bwa-mem2 index -p {params.tmp}/idx "$REF"
                bwa-mem2 mem -t {threads} {params.tmp}/idx {input.r1} \
                  | samtools sort -@ {threads} -T {params.tmp}/srt -o {output.bam} -
            fi
        else
            echo "read_type sem reads para mapear: $RT" >&2; exit 1
        fi
        samtools index -@ {threads} {output.bam}
        rm -rf {params.tmp}
        """


# --------------------------- vRhyme -----------------------------------
rule vrhyme:
    input:
        fna = VIRAL_FNA,
        bam = BAM,
    output:
        bindir = directory(D_BINS + "/vrhyme/{sample}/vRhyme_best_bins_fasta"),
    params:
        outdir = lambda w: f"{D_BINS}/vrhyme/{w.sample}",
    threads: config["threads"]["binning"]
    conda: "../envs/vrhyme.yaml"
    log: D_LOGS + "/03_vrhyme/{sample}.log"
    shell:
        r"""
        set -uo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")"
        echo "=== vRhyme | {wildcards.sample} | $(date) ===" > "$LOG"
        rm -rf {params.outdir}

        # Sem input util -> etapa legitimamente vazia (NAO e erro)
        if [ ! -s {input.fna} ] || [ ! -s {input.bam} ]; then
            mkdir -p {output.bindir}
            echo "[vRhyme] SEM input (fna/bam vazio) -> 0 bins. Etapa pulada legitimamente." >> "$LOG"
            exit 0
        fi

        # Roda vRhyme SEM mascarar o codigo de saida
        vRhyme -i {input.fna} -b {input.bam} -t {threads} -o {params.outdir} >> "$LOG" 2>&1
        rc=$?
        mkdir -p {output.bindir}
        nbins=$(find {output.bindir} -type f \( -name '*.fasta' -o -name '*.fna' \) 2>/dev/null | wc -l)

        # Distingue: sucesso | rodou-mas-0-bins (ok) | erro real
        if [ "$rc" -eq 0 ]; then
            echo "[vRhyme] OK (exit 0) -- rodou ate o fim. Bins formados: $nbins" >> "$LOG"
        elif grep -qiE "no bins|did not (generate|form)|0 bins|not enough|no viral|insufficient" "$LOG"; then
            echo "[vRhyme] OK -- rodou ate o fim, mas NENHUM bin formado (sem sinal de co-cobertura). NAO e erro." >> "$LOG"
        else
            echo "[vRhyme] *** ERRO REAL (exit=$rc) *** -- veja o traceback acima no log. Pipeline ABORTADO de proposito." >> "$LOG"
            exit "$rc"
        fi
        """


# ---------------------------- VAMB ------------------------------------
rule vamb:
    input:
        fna = VIRAL_FNA,
        bam = BAM,
    output:
        bindir = directory(D_BINS + "/vamb/{sample}/bins"),
    params:
        outdir   = lambda w: f"{D_BINS}/vamb/{w.sample}",
        minfasta = config["binning"]["vamb_minfasta"],
    threads: config["threads"]["binning"]
    conda: "../envs/vamb.yaml"
    log: D_LOGS + "/03_vamb/{sample}.log"
    shell:
        r"""
        set -uo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")"
        echo "=== VAMB | {wildcards.sample} | $(date) ===" > "$LOG"
        rm -rf {params.outdir}

        if [ ! -s {input.fna} ] || [ ! -s {input.bam} ]; then
            mkdir -p {output.bindir}
            echo "[VAMB] SEM input (fna/bam vazio) -> 0 bins. Etapa pulada legitimamente." >> "$LOG"
            exit 0
        fi

        vamb bin default --outdir {params.outdir} --fasta {input.fna} \
            --bamfiles {input.bam} --minfasta {params.minfasta} >> "$LOG" 2>&1
        rc=$?
        mkdir -p {output.bindir}
        # VAMB 4.x as vezes escreve em vae_clusters/ em vez de bins/
        if [ ! "$(ls -A {output.bindir} 2>/dev/null)" ] && [ -d {params.outdir}/vae_clusters ]; then
            cp {params.outdir}/vae_clusters/*.fna {output.bindir}/ 2>/dev/null || true
        fi
        nbins=$(find {output.bindir} -type f \( -name '*.fna' -o -name '*.fasta' \) 2>/dev/null | wc -l)

        if [ "$rc" -eq 0 ]; then
            echo "[VAMB] OK (exit 0) -- rodou ate o fim. Bins: $nbins" >> "$LOG"
        elif grep -qiE "no bins|too few|0 bins|no clusters|empty|insufficient" "$LOG"; then
            echo "[VAMB] OK -- rodou ate o fim, mas 0 bins (poucos contigs / sinal fraco). NAO e erro." >> "$LOG"
        else
            echo "[VAMB] *** ERRO REAL (exit=$rc) *** -- veja o log acima. Pipeline ABORTADO de proposito." >> "$LOG"
            exit "$rc"
        fi
        """


rule collect_bins:
    input:
        vr = rules.vrhyme.output.bindir,
        vb = rules.vamb.output.bindir,
    output:
        listing = D_BINS + "/collected/{sample}.bins.tsv",
        outdir  = directory(D_BINS + "/collected/{sample}"),
    params:
        script = f"{config['scripts_dir']}/collect_bins.py",
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        python {params.script} --sample {wildcards.sample} \
            --vrhyme-dir {input.vr} --vamb-dir {input.vb} \
            --out-dir {output.outdir} --listing {output.listing}
        """


# ---------- link dos scaffolds do bin (vRhyme) p/ CheckV no bin ----------
rule link_bins:
    """Une os scaffolds de cada bin com espacador de N's -> 1 seq por bin,
    para o CheckV estimar a qualidade do BIN INTEIRO (recomendacao do vRhyme)."""
    input:
        listing = D_BINS + "/collected/{sample}.bins.tsv",
        bindir  = D_BINS + "/collected/{sample}",
    output:
        fna = LINKED_BINS,
    params:
        script   = f"{config['scripts_dir']}/link_bins.py",
        spacer_n = config.get("binning", {}).get("link_spacer_n", 1000),
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fna})
        python {params.script} --bindir {input.bindir} --listing {input.listing} \
            --spacer-n {params.spacer_n} --out-fna {output.fna}
        """


# ------------------- unidades finais (bin OU contig) ------------------
def filter_inputs(wildcards):
    """Escolhe a fonte por amostra: bins (com reads) ou contigs (sem reads).
    COM reads -> CheckV do BIN unido por N's (CHECKV_BINS_SUM, vMAG inteiro).
    SEM reads -> CheckV por-contig (CHECKV_SUM)."""
    if has_reads(wildcards.sample):
        return {
            "summary": CHECKV_BINS_SUM.format(sample=wildcards.sample),
            "listing": f"{D_BINS}/collected/{wildcards.sample}.bins.tsv",
            "bindir":  f"{D_BINS}/collected/{wildcards.sample}",
        }
    return {
        "summary":   CHECKV_SUM.format(sample=wildcards.sample),
        "viral_fna": VIRAL_FNA.format(sample=wildcards.sample),
    }


def filter_args(wildcards, input):
    if has_reads(wildcards.sample):
        return f"--mode metagenome --listing {input.listing} --bindir {input.bindir}"
    return f"--mode genome --viral-fna {input.viral_fna} --sample {wildcards.sample}"


rule filter_units:
    input:
        unpack(filter_inputs)
    output:
        fna   = PASSED_FNA,
        units = UNITS_TSV,
    params:
        min_vg = config["checkv"]["min_viral_genes"],
        script = f"{config['scripts_dir']}/filter_units.py",
        args   = filter_args,
    conda: "../envs/pyutils.yaml"
    shell:
        r"""
        python {params.script} --summary {input.summary} {params.args} \
            --min-viral-genes {params.min_vg} \
            --out-fna {output.fna} --out-units {output.units}
        """
