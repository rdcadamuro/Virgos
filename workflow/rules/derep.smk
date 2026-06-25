# ======================================================================
# BLOCO 03b — Dereplicação: junta as unidades aprovadas (todas as amostras)
# e agrupa em vOTUs (ANI>=95%, AF>=85%). Fica dentro de 03_bins_vrhyme_vamb/derep.
# ======================================================================

ALL_PASSED_FNA   = expand(PASSED_FNA, sample=SAMPLES)
ALL_PASSED_UNITS = expand(UNITS_TSV, sample=SAMPLES)


rule derep_concat:
    input:
        fnas  = ALL_PASSED_FNA,
        units = ALL_PASSED_UNITS,
    output:
        fna   = D_BINS + "/derep/all_passed.fna",
        units = ALL_UNITS,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fna})
        : > {output.fna}
        for f in {input.fnas}; do [ -s "$f" ] && cat "$f" >> {output.fna} || true; done
        head -n1 {input.units[0]} > {output.units}
        for t in {input.units}; do tail -n+2 "$t" >> {output.units}; done
        """


rule derep_cluster:
    input:
        fna = rules.derep_concat.output.fna,
    output:
        votus    = VOTUS,
        clusters = CLUSTERS,
        votu_map = VOTU_MAP,
    params:
        min_ani = config["derep"]["min_ani"],
        min_af  = config["derep"]["min_af"],
        script  = f"{config['scripts_dir']}/dereplicate.py",
        workdir = D_BINS + "/derep/work",
    threads: 8
    conda: "../envs/derep.yaml"
    log: D_LOGS + "/03_derep.log"
    shell:
        r"""
        set -euo pipefail
        if [ ! -s {input.fna} ]; then
            echo "Nenhuma unidade para dereplicar" > {log}
            : > {output.votus}
            echo -e "representative\tmembers" > {output.clusters}
            echo -e "contig\tvotu_rep" > {output.votu_map}
            exit 0
        fi
        python {params.script} --fna {input.fna} --workdir {params.workdir} \
            --threads {threads} --min-ani {params.min_ani} --min-af {params.min_af} \
            --out-votus {output.votus} --out-clusters {output.clusters} \
            --out-map {output.votu_map} > {log} 2>&1
        """
