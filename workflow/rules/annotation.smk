# ======================================================================
# BLOCO 06/07 — ANOTAÇÃO FUNCIONAL dos vOTUs (sequências QC-aprovadas):
#   RESISTOMA  : ABRicate (multi-DB) -> ARGs / virulência / plasmídeos, dedup por locus
#   MOBILOMA   : consolida vários sinais por vOTU:
#                  - geNomad      : provírus, AMR, conjugação, integrases (já anotado)
#                  - IntegronFinder        : integrons (complete / In0 / CALIN)
#                  - MobileElementFinder   : IS / transposons (MGEdb)
#                  - tRNAscan-SE           : tRNAs (marcam sítios de integração)
#                  - ICEfinder (OPCIONAL)  : ICEs — requer instalação manual
# Roda só sobre os vOTUs. Cada etapa degrada com elegância (saída vazia + nota).
# ======================================================================

_ANN = config.get("annotation", {})
D_ATOOLS = D_ANNOT + "/tools"


# ----------------------------- RESISTOMA ------------------------------
rule abricate_resistome:
    input:
        votus = VOTUS,
    output:
        resistome = D_ANNOT + "/resistome_abricate.tsv",
        by_db     = D_ANNOT + "/resistome_by_db.tsv",
        cp        = D_ANNOT + "/resistome_dedup_checkpoint.json",
    params:
        dbs     = (lambda wildcards: _ANN.get("abricate_dbs", "all")
                   if isinstance(_ANN.get("abricate_dbs", "all"), str)
                   else " ".join(_ANN.get("abricate_dbs"))),
        minid   = _ANN.get("abricate_minid", 80),
        mincov  = _ANN.get("abricate_mincov", 80),
        overlap = _ANN.get("dedup_overlap", 0.5),
        dedup   = f"{config['scripts_dir']}/dedup_abricate.py",
        raw     = D_ANNOT + "/abricate_raw.tsv",
        perdb   = D_ATOOLS + "/abricate",
    threads: 4
    conda: "../envs/abricate.yaml"
    log: D_LOGS + "/06_abricate.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.resistome})" {params.perdb}
        echo "=== ABRicate (resistoma) | $(date) ===" > "$LOG"
        HDR="sequence\tstart\tend\tstrand\tbest_gene\tpct_identity\tpct_coverage\tn_databases\tdatabases\tall_genes\tproduct\tresistance"
        empty() {{
            printf "$HDR\n" > {output.resistome}
            printf "database\tn_raw_hits\tn_loci_supported\n" > {output.by_db}
            printf '{{"n_raw_hits": 0, "n_loci": 0, "verification": {{"residual_overlaps": 0, "status": "PASS"}}, "note": "%s"}}\n' "$1" > {output.cp}
        }}
        if [ ! -s {input.votus} ]; then
            echo "[ABRicate] SEM vOTUs -> tabela vazia. NAO e erro." >> "$LOG"; empty "sem vOTUs"; exit 0
        fi
        if ! command -v abricate >/dev/null 2>&1; then
            echo "[ABRicate] binario 'abricate' nao encontrado -> vazio. Rode: phages.py --block 0" >> "$LOG"; empty "abricate ausente"; exit 0
        fi
        # Resolve a lista de bancos. "all" -> TODOS os bancos que o abricate reporta.
        AVAIL=$(abricate --list 2>/dev/null | awk 'NR>1{{print $1}}' | sort -u)
        DBS="{params.dbs}"
        if [ "$DBS" = "all" ]; then
            DBS=$(echo "$AVAIL" | tr '\n' ' ')
            echo "[ABRicate] modo ALL -> bancos: $DBS" >> "$LOG"
        fi
        RAW="{params.raw}"; : > "$RAW"
        n_ok=0; n_skip=0; RUNDBS=""
        for db in $DBS; do
            if echo "$AVAIL" | grep -qx "$db"; then
                echo "[ABRicate] rodando db=$db (threads={threads})" >> "$LOG"
                # saida por-banco (provenance/comparacao) + concatena no RAW p/ dedup
                if abricate --db "$db" --threads {threads} --minid {params.minid} --mincov {params.mincov} \
                       {input.votus} > "{params.perdb}/$db.tab" 2>> "$LOG"; then
                    cat "{params.perdb}/$db.tab" >> "$RAW"
                    h=$(tail -n +2 "{params.perdb}/$db.tab" 2>/dev/null | wc -l)
                    echo "  -> db=$db: $h hits" >> "$LOG"; n_ok=$((n_ok+1)); RUNDBS="$RUNDBS $db"
                else
                    echo "  (db $db falhou; seguindo)" >> "$LOG"
                fi
            else
                echo "[ABRicate] db ausente: $db (pulado)" >> "$LOG"; n_skip=$((n_skip+1))
            fi
        done
        python {params.dedup} --in "$RAW" --overlap {params.overlap} --dbs-run "$RUNDBS" \
            --out {output.resistome} --out-by-db {output.by_db} \
            --out-checkpoint {output.cp} >> "$LOG" 2>&1
        n=$(tail -n +2 {output.resistome} 2>/dev/null | wc -l)
        st=$(python -c "import json; print(json.load(open('{output.cp}'))['verification']['status'])" 2>/dev/null || echo "?")
        echo "[ABRicate] OK -- rodou ate o fim (exit 0). bancos: $n_ok rodados, $n_skip pulados; loci (dedup): $n; checkpoint: $st" >> "$LOG"
        """


# --------------------- MOBILOMA: geNomad (consolida) ------------------
rule mobilome_genomad:
    input:
        votus = VOTUS,
    output:
        tsv   = D_ATOOLS + "/mobilome_genomad.tsv",
        calls = D_ATOOLS + "/calls_genomad.tsv",
    params:
        script  = f"{config['scripts_dir']}/build_mobilome.py",
        genomad = D_GENOMAD,
        samples = " ".join(SAMPLES) if SAMPLES else "NONE",
        mge     = lambda wildcards: _ANN.get("mge_keywords",
                  "integrase|transposase|recombinase|relaxase|resolvase|mobiliz|conjugat"),
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/06_mobilome_genomad.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS}
        echo "=== Mobiloma/geNomad | $(date) ===" > "$LOG"
        HDR="votu_rep\tsample\ttopology\tis_provirus\tn_amr_genes\tamr_genes\tn_conjugation_genes\tn_integrase\tn_transposase\tn_mge_total"
        CHDR="seq\tstart\tend\telement_type\ttool\tscore"
        if [ ! -s {input.votus} ]; then
            printf "$HDR\n" > {output.tsv}
            printf "$CHDR\n" > {output.calls}
            echo "[mobilome/geNomad] SEM vOTUs -> vazio." >> "$LOG"; exit 0
        fi
        python {params.script} --votus {input.votus} --genomad-dir {params.genomad} \
            --samples {params.samples} --mge-keywords '{params.mge}' \
            --out {output.tsv} --out-calls {output.calls} >> "$LOG" 2>&1
        [ -s {output.calls} ] || printf "$CHDR\n" > {output.calls}
        echo "[mobilome/geNomad] OK -- rodou ate o fim (exit 0)." >> "$LOG"
        """


# --------------------- MOBILOMA: IntegronFinder ----------------------
rule integron_finder:
    input:
        votus = VOTUS,
    output:
        summ = D_ATOOLS + "/integron_summary.tsv",
        done = D_ATOOLS + "/integronfinder/.done",
    params:
        run    = bool(_ANN.get("run_integron_finder", True)),
        outdir = D_ATOOLS + "/integronfinder",
    threads: 4
    conda: "../envs/integronfinder.yaml"
    log: D_LOGS + "/06_integronfinder.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS} {params.outdir}
        echo "=== IntegronFinder | $(date) ===" > "$LOG"
        : > {output.summ}; : > {output.done}
        if [ "{params.run}" != "True" ]; then
            echo "[IntegronFinder] desativado em config (annotation.run_integron_finder=false)." >> "$LOG"; exit 0
        fi
        if [ ! -s {input.votus} ] || ! command -v integron_finder >/dev/null 2>&1; then
            echo "[IntegronFinder] sem vOTUs ou binario ausente -> vazio. NAO e erro." >> "$LOG"; exit 0
        fi
        # limpa apenas os resultados (preserva o sentinel .done que e output)
        find {params.outdir} -mindepth 1 ! -name '.done' -exec rm -rf {{}} + 2>/dev/null || true
        integron_finder --local-max --cpu {threads} --outdir {params.outdir} {input.votus} >> "$LOG" 2>&1 \
            || echo "[IntegronFinder] aviso: execucao retornou erro; seguindo com vazio." >> "$LOG"
        s=$(find {params.outdir} -name '*.summary' 2>/dev/null | head -1 || true)
        [ -n "$s" ] && cp "$s" {output.summ} || true
        n=$(find {params.outdir} -name '*.integrons' 2>/dev/null | wc -l)
        echo "[IntegronFinder] OK -- rodou ate o fim. summary: $([ -s {output.summ} ] && echo presente || echo vazio); .integrons: $n" >> "$LOG"
        """


# ------------------- MOBILOMA: CONJScan / MacSyFinder ----------------
# Alternativa instalavel e host-agnostica ao ICEfinder: detecta sistemas
# conjugativos (T4SS/relaxase) por vOTU. Emite element calls com coordenada.
rule conjscan:
    input:
        votus = VOTUS,
    output:
        calls = D_ATOOLS + "/calls_conjscan.tsv",
    params:
        run    = bool(_ANN.get("run_conjscan", True)),
        script = f"{config['scripts_dir']}/run_conjscan.py",
        models = _ANN.get("conjscan_models", "CONJScan"),
        outdir = D_ATOOLS + "/conjscan",
        overlap = _ANN.get("dedup_overlap", 0.5),
    threads: 4
    conda: "../envs/macsyfinder.yaml"
    log: D_LOGS + "/06_conjscan.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS} {params.outdir}
        echo "=== CONJScan/MacSyFinder | $(date) ===" > "$LOG"
        CHDR="seq\tstart\tend\telement_type\ttool\tscore"
        printf "$CHDR\n" > {output.calls}
        if [ "{params.run}" != "True" ]; then
            echo "[CONJScan] desativado em config (annotation.run_conjscan=false) -> vazio. NAO e erro." >> "$LOG"; exit 0
        fi
        if [ ! -s {input.votus} ] || ! command -v macsyfinder >/dev/null 2>&1; then
            echo "[CONJScan] sem vOTUs ou macsyfinder ausente -> vazio. NAO e erro." >> "$LOG"; exit 0
        fi
        # CONJScan precisa dos modelos instalados localmente. 'macsydata available'
        # lista o INSTALAVEL (remoto), nao o instalado -> checamos o dir de modelos
        # do usuario e instalamos com --user (dir system-wide do env nao e gravavel).
        if [ ! -d "$HOME/.macsyfinder/models/CONJScan" ]; then
            echo "[CONJScan] modelos CONJScan ausentes — instalando (--user)..." >> "$LOG"
            macsydata install --user CONJScan >> "$LOG" 2>&1 \
              || msf_data install --user CONJScan >> "$LOG" 2>&1 \
              || echo "[CONJScan] nao foi possivel instalar os modelos (rede?) -> seguira vazio. NAO e erro." >> "$LOG"
        fi
        python {params.script} --votus {input.votus} --out {output.calls} \
            --models "{params.models}" --workdir {params.outdir} \
            --overlap {params.overlap} --log "$LOG" >> "$LOG" 2>&1 \
            || echo "[CONJScan] aviso: execucao retornou erro; seguindo com vazio." >> "$LOG"
        [ -s {output.calls} ] || printf "$CHDR\n" > {output.calls}
        n=$(tail -n +2 {output.calls} 2>/dev/null | wc -l)
        echo "[CONJScan] OK -- rodou ate o fim (exit 0). sistemas conjugativos (calls): $n" >> "$LOG"
        """


# ----------------- MOBILOMA: MobileElementFinder --------------------
rule mobileelementfinder:
    input:
        votus = VOTUS,
    output:
        csv = D_ATOOLS + "/mefinder.csv",
    params:
        run    = bool(_ANN.get("run_mobileelementfinder", True)),
        outdir = D_ATOOLS + "/mefinder",
    threads: 4
    conda: "../envs/mefinder.yaml"
    log: D_LOGS + "/06_mefinder.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS}
        echo "=== MobileElementFinder | $(date) ===" > "$LOG"
        : > {output.csv}
        if [ "{params.run}" != "True" ]; then
            echo "[mefinder] desativado em config." >> "$LOG"; exit 0
        fi
        if [ ! -s {input.votus} ] || ! command -v mefinder >/dev/null 2>&1; then
            echo "[mefinder] sem vOTUs ou binario ausente -> vazio. NAO e erro." >> "$LOG"; exit 0
        fi
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        ( cd {params.outdir} && mefinder find --contig {input.votus} mefinder ) >> "$LOG" 2>&1 \
            || echo "[mefinder] aviso: execucao retornou erro; seguindo com vazio." >> "$LOG"
        c=$(find {params.outdir} -name 'mefinder*.csv' 2>/dev/null | head -1 || true)
        [ -n "$c" ] && cp "$c" {output.csv} || true
        echo "[mefinder] OK -- rodou ate o fim. csv: $([ -s {output.csv} ] && echo presente || echo vazio)" >> "$LOG"
        """


# ------------------------ MOBILOMA: tRNAscan-SE ----------------------
rule trnascan:
    input:
        votus = VOTUS,
    output:
        tsv = D_ATOOLS + "/trnascan.tsv",
    params:
        run = bool(_ANN.get("run_trnascan", True)),
    threads: 4
    conda: "../envs/trnascan.yaml"
    log: D_LOGS + "/06_trnascan.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS}
        echo "=== tRNAscan-SE | $(date) ===" > "$LOG"
        : > {output.tsv}
        if [ "{params.run}" != "True" ]; then
            echo "[tRNAscan] desativado em config." >> "$LOG"; exit 0
        fi
        if [ ! -s {input.votus} ] || ! command -v tRNAscan-SE >/dev/null 2>&1; then
            echo "[tRNAscan] sem vOTUs ou binario ausente -> vazio. NAO e erro." >> "$LOG"; exit 0
        fi
        # tRNAscan-SE ABORTA se o arquivo -o ja existir -> escrevemos num path
        # temporario inexistente e depois copiamos. -G = modo geral (qualquer dominio).
        TMP="{output.tsv}.raw"; rm -f "$TMP"
        tRNAscan-SE -G --thread {threads} -q -o "$TMP" {input.votus} >> "$LOG" 2>&1 \
            || echo "[tRNAscan] aviso: execucao retornou erro; seguindo com vazio." >> "$LOG"
        if [ -s "$TMP" ]; then cp "$TMP" {output.tsv}; fi
        rm -f "$TMP"
        n=$(grep -cve '^Sequence\|^Name\|^---\|^$' {output.tsv} 2>/dev/null || true)
        echo "[tRNAscan] OK -- rodou ate o fim. linhas de tRNA: ${{n:-0}}" >> "$LOG"
        """


# ------------------- MOBILOMA: ICEfinder (OPCIONAL) -----------------
rule icefinder:
    input:
        votus = VOTUS,
    output:
        tsv = D_ATOOLS + "/icefinder.tsv",
    params:
        run  = bool(_ANN.get("run_icefinder", False)),
        path = _ANN.get("icefinder_path", ""),
    log: D_LOGS + "/06_icefinder.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {D_ATOOLS}
        echo "=== ICEfinder (opcional) | $(date) ===" > "$LOG"
        : > {output.tsv}
        if [ "{params.run}" != "True" ] || [ -z "{params.path}" ]; then
            echo "[ICEfinder] desativado/ausente (nao e conda; instale manualmente e aponte annotation.icefinder_path). Pulado." >> "$LOG"
            exit 0
        fi
        echo "[ICEfinder] usando: {params.path} (integracao depende da sua instalacao manual)" >> "$LOG"
        # Gancho: a saida do ICEfinder varia por instalacao; consolidate_mobilome.py
        # le {output.tsv} se voce gravar aqui o resultado no formato esperado.
        """


# ----------------------- MOBILOMA: consolidação ----------------------
rule mobilome_consolidate:
    input:
        votus     = VOTUS,
        genomad   = rules.mobilome_genomad.output.tsv,
        gcalls    = rules.mobilome_genomad.output.calls,
        integron  = rules.integron_finder.output.summ,
        intdone   = rules.integron_finder.output.done,
        mefinder  = rules.mobileelementfinder.output.csv,
        conjscan  = rules.conjscan.output.calls,
        trnascan  = rules.trnascan.output.tsv,
        icefinder = rules.icefinder.output.tsv,
    output:
        mob     = D_ANNOT + "/mobilome_table.tsv",
        regions = D_ANNOT + "/mobile_element_regions.tsv",
        cp      = D_ANNOT + "/mobilome_dedup_checkpoint.json",
    params:
        script  = f"{config['scripts_dir']}/consolidate_mobilome.py",
        intdir  = D_ATOOLS + "/integronfinder",
        overlap = _ANN.get("dedup_overlap", 0.5),
    conda: "../envs/pyutils.yaml"
    log: D_LOGS + "/06_mobilome.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" "$(dirname {output.mob})"
        echo "=== Mobiloma (consolidacao + dedup 2 camadas) | $(date) ===" > "$LOG"
        python {params.script} \
            --votus {input.votus} --genomad {input.genomad} \
            --integron {input.integron} --integron-dir {params.intdir} \
            --mefinder {input.mefinder} \
            --trnascan {input.trnascan} --icefinder {input.icefinder} \
            --genomad-calls {input.gcalls} --conjscan-calls {input.conjscan} \
            --icefinder-calls {input.icefinder} \
            --overlap {params.overlap} \
            --out {output.mob} --out-regions {output.regions} \
            --out-checkpoint {output.cp} >> "$LOG" 2>&1
        n=$(tail -n +2 {output.mob} 2>/dev/null | wc -l)
        r=$(tail -n +2 {output.regions} 2>/dev/null | wc -l)
        st=$(python -c "import json,sys; print(json.load(open('{output.cp}'))['verification']['status'])" 2>/dev/null || echo "?")
        echo "[mobilome] OK -- rodou ate o fim (exit 0). vOTUs anotados: $n; regioes MGE nao-redundantes: $r; dedup checkpoint: $st" >> "$LOG"
        """
