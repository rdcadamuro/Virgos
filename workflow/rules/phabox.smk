# ======================================================================
# BLOCO 04 — PhaBox2 sobre os vOTUs. Módulo PhaTYP -> lifestyle.
# ======================================================================

rule phabox:
    input:
        votus = VOTUS,
        db    = config["databases"]["phabox"],
    output:
        phatyp = D_PHABOX + "/final_prediction/phatyp_prediction.tsv",
    params:
        outdir = D_PHABOX,
        task   = config["phabox"]["task"],
        patch  = f"{config['scripts_dir']}/patch_phabox2.py",
    threads: config["threads"]["phabox"]
    conda: "../envs/phabox.yaml"
    log: D_LOGS + "/04_phabox.log"
    shell:
        r"""
        set -euo pipefail
        LOG="{log}"; mkdir -p "$(dirname "$LOG")" {params.outdir}
        echo "=== PhaBox2 | $(date) ===" > "$LOG"
        if [ ! -s {input.votus} ]; then
            echo "[PhaBox2] SEM vOTUs de entrada -> tabela vazia. NAO e erro." >> "$LOG"
            mkdir -p $(dirname {output.phatyp})
            echo -e "Accession\tTYPE\tScore" > {output.phatyp}
            exit 0
        fi
        # Garante compatibilidade pandas>=2 (idempotente; ver patch_phabox2.py)
        python {params.patch} >> "$LOG" 2>&1 || true
        # AJUSTAR conforme versao do PhaBox2 (--dbdir / --outpth podem variar)
        phabox2 --task {params.task} --contigs {input.votus} \
            --outpth {params.outdir} --dbdir {input.db} --threads {threads} \
            >> "$LOG" 2>&1
        if [ ! -s {output.phatyp} ]; then
            alt=$(find {params.outdir} -name 'phatyp_prediction.*' | head -n1 || true)
            [ -n "$alt" ] && cp "$alt" {output.phatyp} || true
        fi
        # --- VERIFICACAO: chegou aqui => PhaBox2 COMPLETOU (exit 0) ---
        n=$(tail -n +2 {output.phatyp} 2>/dev/null | wc -l)
        echo "[PhaBox2] OK -- rodou ate o fim (exit 0). Predicoes de lifestyle: $n" >> "$LOG"
        """
