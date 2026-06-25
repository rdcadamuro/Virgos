# ======================================================================
# BLOCO 0 (setup) — envs isolados por ferramenta + download dos bancos.
#
#   snakemake --use-conda --conda-create-envs-only setup_envs   # constrói envs
#   snakemake --use-conda setup_dbs                             # baixa DBs
#
# Cada ferramenta vive em SEU env (workflow/envs/*.yaml) -> sem conflito.
# O download de cada DB roda DENTRO do env da própria ferramenta.
# ======================================================================

ENVS = ["genomad", "checkv", "vrhyme", "vamb", "mapping",
        "derep", "phabox", "vitap", "abundance", "abricate",
        "integronfinder", "mefinder", "trnascan", "macsyfinder", "pyutils"]


rule env_build_marker:
    """Job trivial só para forçar a criação do env (com --conda-create-envs-only,
    o env é construído mas o shell não roda)."""
    output: touch(SETUP_DBDIR + "/.envs/{env}.built")
    conda: "../envs/{env}.yaml"
    shell: "echo env pronto: {wildcards.env}"


rule setup_envs:
    input: expand(SETUP_DBDIR + "/.envs/{env}.built", env=ENVS)


# ----------------------------- bancos --------------------------------
rule setup_genomad_db:
    output: touch(SETUP_DBDIR + "/.genomad_db.done")
    conda: "../envs/genomad.yaml"
    log: D_LOGS + "/setup/genomad_db.log"
    shell:
        r"""
        mkdir -p {SETUP_DBDIR}
        if [ -d {SETUP_DBDIR}/genomad_db ]; then
            echo "geNomad DB já existe" | tee {log}
        else
            genomad download-database {SETUP_DBDIR} > {log} 2>&1
        fi
        echo ">> geNomad DB: {SETUP_DBDIR}/genomad_db"
        """


rule setup_checkv_db:
    output: touch(SETUP_DBDIR + "/.checkv_db.done")
    conda: "../envs/checkv.yaml"
    log: D_LOGS + "/setup/checkv_db.log"
    shell:
        r"""
        mkdir -p {SETUP_DBDIR}
        if ls -d {SETUP_DBDIR}/checkv-db* >/dev/null 2>&1; then
            echo "CheckV DB já existe" | tee {log}
        else
            checkv download_database {SETUP_DBDIR} > {log} 2>&1
        fi
        echo ">> CheckV DB: $(ls -d {SETUP_DBDIR}/checkv-db* 2>/dev/null | head -1)"
        """


rule setup_phabox_db:
    """PhaBox2 DB: baixa automaticamente do GitHub Releases (v2.2, ~3 GB zip)."""
    output: touch(SETUP_DBDIR + "/.phabox_db.done")
    log: D_LOGS + "/setup/phabox_db.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SETUP_DBDIR}
        DEST="{SETUP_DBDIR}/phabox_db_v2"

        if [ -d "$DEST" ] && [ "$(ls -A "$DEST" 2>/dev/null | wc -l)" -gt 5 ]; then
            echo "PhaBox2 DB já existe em $DEST" | tee {log}
            exit 0
        fi

        # URL oficial do GitHub Releases (PhaBOX v2.2 — ICTV 2024)
        DL_URL="https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_2.zip"
        TMP_ZIP="{SETUP_DBDIR}/phabox_db_v2_2.zip"

        # Aproveita zip já existente (e.g., de run anterior interrompido)
        if file "$TMP_ZIP" 2>/dev/null | grep -qiE "zip|compress"; then
            echo "==> Zip já existe e é válido, pulando download." | tee -a {log}
        else
            echo "==> Baixando PhaBox2 DB (~3 GB, pode demorar)..." | tee -a {log}
            curl -L --progress-bar -o "$TMP_ZIP" "$DL_URL" 2>> {log} || \
                wget -L -q --show-progress -O "$TMP_ZIP" "$DL_URL" >> {log} 2>&1
        fi

        # Valida que é zip real (não HTML de erro)
        if ! file "$TMP_ZIP" | grep -qiE "zip|compress"; then
            echo "ERRO: arquivo baixado não é zip. Primeiras linhas:" | tee -a {log}
            head -3 "$TMP_ZIP" | tee -a {log}
            rm -f "$TMP_ZIP"
            exit 1
        fi

        echo "==> Extraindo..." | tee -a {log}
        python3 -m zipfile -e "$TMP_ZIP" "{SETUP_DBDIR}" >> {log} 2>&1
        rm -f "$TMP_ZIP"

        # Garante que o diretório final se chama phabox_db_v2
        if [ ! -d "$DEST" ]; then
            extracted=$(find "{SETUP_DBDIR}" -maxdepth 1 -type d -name "phabox_db*" | head -1)
            [ -n "$extracted" ] && mv "$extracted" "$DEST"
        fi

        echo ">> PhaBox2 DB: $DEST" | tee -a {log}
        """


rule setup_vitap_db:
    """VITAP DB: baixa pre-built DB do Figshare (25426159) via API."""
    output: touch(SETUP_DBDIR + "/.vitap_db.done")
    log: D_LOGS + "/setup/vitap_db.log"
    params:
        db_dir = SETUP_DBDIR + "/ViTAP_db",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SETUP_DBDIR}

        if [ -d "{params.db_dir}" ] && [ "$(ls -A "{params.db_dir}" 2>/dev/null | wc -l)" -gt 3 ]; then
            echo "VITAP DB já existe em {params.db_dir}" | tee {log}
            exit 0
        fi

        # Consulta Figshare API para pegar URL de download direta
        FIGSHARE_ID="25426159"
        echo "==> Consultando Figshare para VITAP DB (artigo $FIGSHARE_ID)..." | tee {log}

        DL_URL=$(python3 - <<'PYEOF'
import urllib.request, json, sys
try:
    url = "https://api.figshare.com/v2/articles/25426159/files"
    data = json.loads(urllib.request.urlopen(url).read())
    for f in data:
        name = f.get("name", "").lower()
        if name.endswith(".zip") or "db" in name or "vitap" in name:
            print(f["download_url"])
            sys.exit(0)
    if data:
        print(data[0]["download_url"])
except Exception as e:
    pass
PYEOF
)

        if [ -z "$DL_URL" ]; then
            echo "AVISO: Figshare API nao respondeu. Tente executar block 0 novamente." | tee -a {log}
            exit 1
        fi

        echo "URL obtida: $DL_URL" | tee -a {log}
        TMP_ZIP="{SETUP_DBDIR}/vitap_db.zip"

        echo "==> Baixando VITAP DB (~620 MB)..." | tee -a {log}
        curl -L --progress-bar -o "$TMP_ZIP" "$DL_URL" 2>> {log} || \
            wget -L -q --show-progress -O "$TMP_ZIP" "$DL_URL" >> {log} 2>&1

        if ! file "$TMP_ZIP" | grep -qiE "zip|compress"; then
            echo "ERRO: arquivo baixado nao e zip." | tee -a {log}
            head -3 "$TMP_ZIP" | tee -a {log}
            rm -f "$TMP_ZIP"; exit 1
        fi

        echo "==> Extraindo VITAP DB..." | tee -a {log}
        mkdir -p "{params.db_dir}"
        python3 -m zipfile -e "$TMP_ZIP" "{params.db_dir}" >> {log} 2>&1
        rm -f "$TMP_ZIP"

        echo ">> VITAP DB: {params.db_dir}" | tee -a {log}
        """


rule setup_abricate_db:
    """ABRicate ja traz TODOS os bancos embutidos e indexados (nao ha download).
    Aqui VALIDAMOS (abricate --check) e gravamos um MANIFESTO (abricate --list)
    p/ proveniencia/reprodutibilidade: qual versao/data de cada banco foi usada.
    Para ATUALIZAR um banco a versao mais recente da fonte (opcional, nao-
    deterministico): 'abricate-get_db --db <nome> --force' dentro do env abricate."""
    output:
        done     = touch(SETUP_DBDIR + "/.abricate_db.done"),
        manifest = SETUP_DBDIR + "/abricate_db_manifest.tsv",
    conda: "../envs/abricate.yaml"
    log: D_LOGS + "/setup/abricate_db.log"
    shell:
        r"""
        set -uo pipefail
        mkdir -p {SETUP_DBDIR} "$(dirname {log})"
        echo "=== ABRicate DB setup | $(date) ===" > {log}
        if ! command -v abricate >/dev/null 2>&1; then
            echo "[abricate] binario ausente -> manifesto vazio (rode block 0)." | tee -a {log}
            printf "DATABASE\tSEQUENCES\tDBTYPE\tDATE\n" > {output.manifest}; exit 0
        fi
        abricate --check >> {log} 2>&1 || echo "[abricate] --check reportou aviso (ver log)" >> {log}
        abricate --list > {output.manifest} 2>> {log} || printf "DATABASE\tSEQUENCES\tDBTYPE\tDATE\n" > {output.manifest}
        n=$(tail -n +2 {output.manifest} 2>/dev/null | wc -l)
        echo ">> ABRicate: $n bancos disponiveis (embutidos) -> {output.manifest}" | tee -a {log}
        """


rule setup_dbs:
    input:
        SETUP_DBDIR + "/.genomad_db.done",
        SETUP_DBDIR + "/.checkv_db.done",
        SETUP_DBDIR + "/.phabox_db.done",
        SETUP_DBDIR + "/.vitap_db.done",
        SETUP_DBDIR + "/.abricate_db.done",
    run:
        import os
        vitap_tools = os.path.join(SETUP_DBDIR, "ViTAP_tools")
        print("\n==================== BANCOS DE DADOS ====================")
        print(f"geNomad : {SETUP_DBDIR}/genomad_db")
        print(f"CheckV  : {SETUP_DBDIR}/checkv-db-*")
        print(f"PhaBox2 : {SETUP_DBDIR}/phabox_db_v2")
        print(f"ViTAP   : {SETUP_DBDIR}/ViTAP_db  (script: {vitap_tools}/ViTAP.py)")
        print(f"ABRicate: embutidos -> manifesto {SETUP_DBDIR}/abricate_db_manifest.tsv")
        print("\nCaminhos já atualizados em config/config.yaml -> databases")
        print("=========================================================")
