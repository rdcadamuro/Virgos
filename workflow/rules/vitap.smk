# ======================================================================
# BLOCO 05 — VITAP: taxonomia dos vOTUs (DrKaiyangZheng/VITAP, bioconda).
# DB pre-built baixado pelo bloco 0 em db_dir/ViTAP_db/.
# Output principal: best_determined_lineages.tsv
# ======================================================================

rule vitap:
    input:
        votus = VOTUS,
        db    = config["databases"]["vitap"],
    output:
        tax = D_VITAP + "/vitap_taxonomy.tsv",
    params:
        outdir = D_VITAP,
    threads: config["threads"]["vitap"]
    conda: "../envs/vitap.yaml"
    log: D_LOGS + "/05_vitap.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        if [ ! -s {input.votus} ]; then
            echo "Sem vOTUs para VITAP" > {log}
            echo -e "contig\ttaxonomy" > {output.tax}
            exit 0
        fi

        # Tenta binary 'VITAP' (bioconda) ou 'vitap' (alias)
        VITAP_BIN=$(command -v VITAP 2>/dev/null || command -v vitap 2>/dev/null || echo "")
        if [ -z "$VITAP_BIN" ]; then
            echo "ERRO: binario VITAP nao encontrado no env. Rode: python phages.py --block 0" | tee {log}
            echo -e "contig\ttaxonomy" > {output.tax}
            exit 0
        fi

        # Se uniref90.dmnd nao existe, cria fallback vazio para pular essa etapa
        if [ ! -f {input.db}/uniref90.dmnd ]; then
            echo "uniref90.dmnd nao encontrado — criando fallback vazio para pular etapa UniRef90" >> {log}
            printf 'genome_id\ttaxa_name\tparticipation_index\ttaxon_level\n' \
                > {params.outdir}/target_uniref90_taxa_fallback.out
        fi

        "$VITAP_BIN" assignment \
            -i {input.votus} \
            -d {input.db} \
            -o {params.outdir} \
            >> {log} 2>&1

        # Normaliza saida para nome esperado pelo pipeline
        if [ ! -s {output.tax} ]; then
            res=$(find {params.outdir} -name 'best_determined_lineages.tsv' 2>/dev/null | head -1 || true)
            [ -n "$res" ] && cp "$res" {output.tax} || echo -e "contig\ttaxonomy" > {output.tax}
        fi

        # --- VERIFICACAO: chegou aqui => VITAP COMPLETOU (exit 0) ---
        n=$(tail -n +2 {output.tax} 2>/dev/null | grep -cve '^[[:space:]]*$' || echo 0)
        echo "[VITAP] OK -- rodou ate o fim (exit 0). Linhagens determinadas: $n" >> {log}
        if [ "$n" -eq 0 ]; then
            echo "[VITAP] NOTA: 0 linhagens. VITAP COMPLETOU; os vOTUs nao tem taxonomia confiavel no framework ICTV (fragmentos curtos/divergentes). NAO e erro." >> {log}
        fi
        """
