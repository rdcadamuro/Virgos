# ======================================================================
# ENRIQUECIMENTO — ARG/MGE estão sobre-representados na fração VIRAL
# (prófagos/vOTUs) vs. o GENOMA COMPLETO (background)?
#
#   build_gene_table : rotula TODOS os genes do geNomad (universo) como
#                      viral / ARG / MGE (background "de graça" do geNomad,
#                      + augmentações opcionais por coordenada).
#   enrichment_test  : permutação estratificada por amostra (+ hipergeométrica)
#                      -> enrichment_results.tsv + checkpoint (guardrails).
# Requer os intermediários do geNomad (annotate/_genes.tsv, virus_summary,
# provirus). Se ausentes (ex.: --cleanup agressivo) -> degrada com nota.
# ======================================================================

_ANN_E = config.get("annotation", {})
GENE_TABLE = D_ANNOT + "/enrichment/genes_{sample}.tsv"
INTEGRON_HITS = D_ANNOT + "/enrichment/integron_hits_{sample}.tsv"
IS_HITS       = D_ANNOT + "/enrichment/is_hits_{sample}.tsv"
CONJ_HITS     = D_ANNOT + "/enrichment/conjscan_hits_{sample}.tsv"


# IntegronFinder no GENOMA INTEIRO -> hits de integron por coordenada, para
# rotular genes como is_integron no teste de enriquecimento (categoria que o
# geNomad nao anota). Roda 1x por amostra no assembly completo.
rule integron_genome:
    input:
        assembly = s_assembly,
    output:
        hits = INTEGRON_HITS,
    params:
        run    = bool(_ANN_E.get("run_integron_genome", True)),
        outdir = D_ANNOT + "/enrichment/integronfinder_genome/{sample}",
        parser = f"{config['scripts_dir']}/integron_to_hits.py",
    threads: 4
    conda: "../envs/integronfinder.yaml"
    log: D_LOGS + "/06_integron_genome/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.hits})"
        echo "=== IntegronFinder (genoma inteiro) | {wildcards.sample} | $(date) ===" > "$LOG"
        printf "contig\tstart\tend\n" > {output.hits}
        if [ "{params.run}" != "True" ] || [ ! -s "{input.assembly}" ] || ! command -v integron_finder >/dev/null 2>&1; then
            echo "[integron_genome] desativado / sem assembly / binario ausente -> sem hits. NAO e erro." >> "$LOG"; exit 0
        fi
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        integron_finder --cpu {threads} --outdir {params.outdir} "{input.assembly}" >> "$LOG" 2>&1 \
            || echo "[integron_genome] aviso: IntegronFinder retornou erro; seguindo sem hits." >> "$LOG"
        D=$(find {params.outdir} -maxdepth 1 -type d -name 'Results_Integron_Finder_*' | head -1 || true)
        [ -n "$D" ] && python {params.parser} --integron-dir "$D" --out {output.hits} >> "$LOG" 2>&1 || true
        n=$(tail -n +2 {output.hits} 2>/dev/null | wc -l)
        echo "[integron_genome] OK -- rodou ate o fim. integrons no genoma: $n" >> "$LOG"
        """


# MobileElementFinder no GENOMA INTEIRO -> hits de IS/transposon por coordenada,
# para rotular genes como is_IS no enriquecimento (categoria que o geNomad nao
# anota). Roda 1x por amostra no assembly completo (1 chamada mefinder).
rule mefinder_genome:
    input:
        assembly = s_assembly,
    output:
        hits = IS_HITS,
    params:
        run    = bool(_ANN_E.get("run_is_genome", True)),
        outdir = D_ANNOT + "/enrichment/mefinder_genome/{sample}",
        parser = f"{config['scripts_dir']}/mefinder_to_hits.py",
    threads: 4
    conda: "../envs/mefinder.yaml"
    log: D_LOGS + "/06_is_genome/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.hits})"
        echo "=== MobileElementFinder (genoma inteiro) | {wildcards.sample} | $(date) ===" > "$LOG"
        printf "contig\tstart\tend\n" > {output.hits}
        if [ "{params.run}" != "True" ] || [ ! -s "{input.assembly}" ] || ! command -v mefinder >/dev/null 2>&1; then
            echo "[is_genome] desativado / sem assembly / mefinder ausente -> sem hits. NAO e erro." >> "$LOG"; exit 0
        fi
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        AB="$(readlink -f "{input.assembly}")"
        ( cd {params.outdir} && mefinder find --contig "$AB" mefinder ) >> "$LOG" 2>&1 \
            || echo "[is_genome] aviso: mefinder retornou erro; seguindo sem hits." >> "$LOG"
        C=$(find {params.outdir} -name 'mefinder*.csv' 2>/dev/null | head -1 || true)
        [ -n "$C" ] && python {params.parser} --csv "$C" --out {output.hits} >> "$LOG" 2>&1 || true
        n=$(tail -n +2 {output.hits} 2>/dev/null | wc -l)
        echo "[is_genome] OK -- rodou ate o fim. IS/transposons no genoma: ${{n:-0}}" >> "$LOG"
        """


# CONJScan/MacSyFinder no GENOMA INTEIRO -> hits dos genes de sistemas
# conjugativos, para rotular genes como is_conjugative_system no enriquecimento.
# EFICIENTE: 1 prodigal + 1 macsyfinder (db-type unordered) no genoma todo.
rule conjscan_genome:
    input:
        assembly = s_assembly,
    output:
        hits = CONJ_HITS,
    params:
        run    = bool(_ANN_E.get("run_conjscan_genome", True)),
        script = f"{config['scripts_dir']}/run_conjscan_genome.py",
        models = _ANN_E.get("conjscan_models", "CONJScan"),
        outdir = D_ANNOT + "/enrichment/conjscan_genome/{sample}",
    threads: 4
    conda: "../envs/macsyfinder.yaml"
    log: D_LOGS + "/06_conjscan_genome/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.hits})"
        echo "=== CONJScan (genoma inteiro) | {wildcards.sample} | $(date) ===" > "$LOG"
        printf "contig\tstart\tend\n" > {output.hits}
        if [ "{params.run}" != "True" ] || [ ! -s "{input.assembly}" ] || ! command -v macsyfinder >/dev/null 2>&1; then
            echo "[conjscan_genome] desativado / sem assembly / macsyfinder ausente -> sem hits. NAO e erro." >> "$LOG"; exit 0
        fi
        # CONJScan precisa dos modelos instalados localmente. 'macsydata available'
        # lista o que e INSTALAVEL (remoto), nao o instalado -> checamos o diretorio
        # de modelos do usuario (~/.macsyfinder/models) e instalamos com --user
        # (o dir system-wide do env nao e gravavel).
        if [ ! -d "$HOME/.macsyfinder/models/CONJScan" ]; then
            echo "[conjscan_genome] modelos CONJScan ausentes — instalando (--user)..." >> "$LOG"
            macsydata install --user CONJScan >> "$LOG" 2>&1 \
              || msf_data install --user CONJScan >> "$LOG" 2>&1 \
              || echo "[conjscan_genome] nao foi possivel instalar (rede?) -> seguira vazio. NAO e erro." >> "$LOG"
        fi
        python {params.script} --assembly "{input.assembly}" --out {output.hits} \
            --models "{params.models}" --workdir {params.outdir} --log "$LOG" >> "$LOG" 2>&1 \
            || echo "[conjscan_genome] aviso: execucao retornou erro; seguindo sem hits." >> "$LOG"
        [ -s {output.hits} ] || printf "contig\tstart\tend\n" > {output.hits}
        n=$(tail -n +2 {output.hits} 2>/dev/null | wc -l)
        echo "[conjscan_genome] OK -- rodou ate o fim. genes de sistema conjugativo: ${{n:-0}}" >> "$LOG"
        """


rule build_gene_table:
    input:
        viral    = VIRAL_FNA,          # garante que o geNomad rodou para a amostra
        integron = INTEGRON_HITS,      # is_integron (IntegronFinder no genoma inteiro)
        is_hits  = IS_HITS,            # is_IS (MobileElementFinder no genoma inteiro)
        conj     = CONJ_HITS,          # is_conjugative_system (CONJScan no genoma inteiro)
        phatyp   = D_PHABOX + "/final_prediction/phatyp_prediction.tsv",  # lifestyle (temperate/virulent)
    output:
        tsv = GENE_TABLE,
    params:
        script = f"{config['scripts_dir']}/build_gene_table.py",
        gdir   = D_GENOMAD,
        mge    = lambda wildcards: _ANN_E.get("mge_keywords",
                 "integrase|transposase|recombinase|relaxase|resolvase|mobiliz|conjugat|insertion sequence"),
        vset    = _ANN_E.get("enrichment_viral_set", "temperate_phage"),
        lenient = str(_ANN_E.get("phage_filter_lenient", True)).lower(),
        phage   = lambda wildcards: _ANN_E.get("phage_taxa_classes", "Caudoviricetes"),
        euk     = lambda wildcards: _ANN_E.get("eukaryotic_taxa",
                  "Herviviricetes|Megaviricetes|Revtraviricetes|Artverviricota|Negarnaviricota"),
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/06_gene_table/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.tsv})"
        echo "=== gene table | {wildcards.sample} | $(date) ===" > "$LOG"
        B="{params.gdir}/{wildcards.sample}"
        GENES="$B/{wildcards.sample}_annotate/{wildcards.sample}_genes.tsv"
        VS="$B/{wildcards.sample}_summary/{wildcards.sample}_virus_summary.tsv"
        PV="$B/{wildcards.sample}_find_proviruses/{wildcards.sample}_provirus.tsv"
        if [ ! -s "$GENES" ]; then
            echo "[gene_table] geNomad annotate genes ausente ($GENES) — provavelmente removido por --cleanup. Tabela vazia (enriquecimento sem dados)." >> "$LOG"
            printf "gene\tsample\tcontig\tstart\tend\tis_viral\tis_arg\tis_conjugation\tis_integrase\tis_transposase\tis_recombinase\tis_integron\tis_conjugative_system\tis_IS\tis_mge\n" > {output.tsv}
            exit 0
        fi
        python {params.script} --genes "$GENES" --virus-summary "$VS" --provirus "$PV" \
            --sample {wildcards.sample} --mge-keywords '{params.mge}' \
            --integron-hits {input.integron} \
            --is-hits {input.is_hits} \
            --conjscan-hits {input.conj} \
            --phatyp {input.phatyp} \
            --viral-set {params.vset} --phage-lenient {params.lenient} \
            --phage-classes '{params.phage}' --eukaryotic-classes '{params.euk}' \
            --out {output.tsv} >> "$LOG" 2>&1
        echo "[gene_table] OK -- rodou ate o fim (exit 0)." >> "$LOG"
        """


rule enrichment_test:
    input:
        tables = expand(GENE_TABLE, sample=SAMPLES),
    output:
        tsv = D_ANNOT + "/enrichment_results.tsv",
        cp  = D_ANNOT + "/enrichment_checkpoint.json",
    params:
        script = f"{config['scripts_dir']}/permutation_enrichment.py",
        merged = D_ANNOT + "/enrichment/all_genes.tsv",
        perms  = _ANN_E.get("enrichment_permutations", 10000),
        minev  = _ANN_E.get("enrichment_min_events", 5),
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/06_enrichment.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.tsv})"
        echo "=== enriquecimento (permutacao) | $(date) ===" > "$LOG"
        M="{params.merged}"; : > "$M"; first=1
        for t in {input.tables}; do
            [ -s "$t" ] || continue
            if [ "$first" = "1" ]; then cat "$t" > "$M"; first=0; else tail -n +2 "$t" >> "$M"; fi
        done
        HDR="category\tstatus\tN_genes\tK_positive\tn_viral_genes\tn_nonviral_genes\tk_observed_viral_positive\texpected_viral_positive\tenrichment_ratio\tp_enrichment_perm\tp_depletion_perm\tp_enrichment_hypergeom\tpermutations"
        if [ ! -s "$M" ] || [ "$(tail -n +2 "$M" 2>/dev/null | wc -l)" = "0" ]; then
            echo "[enrichment] sem genes (geNomad annotate ausente) -> sem teste. NAO e erro." >> "$LOG"
            printf "$HDR\n" > {output.tsv}
            printf '{{"results": [], "checkpoint": {{"any_underpowered": false, "any_no_background": true, "note": "sem genes"}}}}\n' > {output.cp}
            exit 0
        fi
        # sem --categories: o motor auto-detecta TODAS as colunas is_* (exceto is_viral)
        python {params.script} --genes "$M" \
            --permutations {params.perms} --min-events {params.minev} \
            --out-tsv {output.tsv} --out-json {output.cp} >> "$LOG" 2>&1
        echo "[enrichment] OK -- rodou ate o fim (exit 0)." >> "$LOG"
        """


rule annotation_report:
    """Relatorio consolidado: TSV geral por vOTU + LOG estatistico (enriquecimento
    por categoria + dedup + totais). Roda apos mobiloma + resistoma + enriquecimento."""
    input:
        mob = D_ANNOT + "/mobilome_table.tsv",
        res = D_ANNOT + "/resistome_abricate.tsv",
        reg = D_ANNOT + "/mobile_element_regions.tsv",
        enr = D_ANNOT + "/enrichment_results.tsv",
    output:
        tsv = D_ANNOT + "/annotation_master.tsv",
        rep = D_ANNOT + "/annotation_stats.log",
    params:
        script = f"{config['scripts_dir']}/annotation_report.py",
        annot  = D_ANNOT,
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/06_annotation_report.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} --annot-dir {params.annot} \
            --out-tsv {output.tsv} --out-log {output.rep} > "{log}" 2>&1
        echo "[annotation_report] OK -- rodou ate o fim (exit 0)." >> "{log}"
        """
