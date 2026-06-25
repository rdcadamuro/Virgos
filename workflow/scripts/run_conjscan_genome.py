#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONJScan/MacSyFinder no GENOMA INTEIRO -> TSV de HITS por coordenada
(contig, start, end) dos genes que compoem sistemas conjugativos, usado por
build_gene_table.py (--conjscan-hits) para rotular genes como is_conjugative_system
no teste de enriquecimento.

Sem isto, is_conjugative_system fica estruturalmente vazio (0 -> UNDERPOWERED):
o CONJScan so rodava por contig nos vOTUs (sem background nao-viral, sem
enriquecimento).

EFICIENTE: 1 prodigal + 1 macsyfinder no genoma todo, com --db-type unordered
(modo recomendado p/ assemblies multi-contig/incompletos), em vez de varrer
centenas de contigs um a um.

BULLETPROOF: sem prodigal/macsyfinder, sem modelos CONJScan, sem proteinas, ou
erro do macsyfinder -> grava TSV so com header + nota; NUNCA falha dura (exit 0).
"""
import argparse
import os
import shutil
import subprocess
import sys


def log(msg, fh):
    fh.write(msg + "\n")
    fh.flush()


def run_prodigal_all(assembly, faa, gff, logfh):
    """prodigal no genoma inteiro (-p meta, robusto p/ draft). True se gerou .faa."""
    try:
        subprocess.run(["prodigal", "-i", assembly, "-a", faa, "-f", "gff",
                        "-o", gff, "-p", "meta", "-q"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        log(f"  [prodigal] falhou: {e}", logfh)
        return False
    return os.path.exists(faa) and os.path.getsize(faa) > 0


def faa_coords(faa):
    """prodigal grava o header como '>{protid} # start # end # strand # ...'.
    -> {protid: (contig, start, end)}; contig = protid sem o ultimo _<idx>."""
    coords = {}
    with open(faa, errors="ignore") as fh:
        for ln in fh:
            if not ln.startswith(">"):
                continue
            parts = ln[1:].split("#")
            pid = parts[0].strip()
            if len(parts) < 3 or not pid:
                continue
            try:
                s, e = int(parts[1].strip()), int(parts[2].strip())
            except ValueError:
                continue
            contig = pid.rsplit("_", 1)[0]
            coords[pid] = (contig, min(s, e), max(s, e))
    return coords


def parse_macsy_hits(out_dir):
    """Le o resultado do macsyfinder -> set de hit_id (id da proteina prodigal).
    Em --db-type unordered o arquivo e all_systems.tsv, que REPETE o cabecalho
    (linha que comeca com 'replicon') a cada sistema -> tratamos cada uma como
    novo cabecalho (reseta o indice de colunas)."""
    hits = set()
    for fn in ("all_systems.tsv", "best_solution.tsv", "all_best_solutions.tsv"):
        path = os.path.join(out_dir, fn)
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            continue
        idx = None
        with open(path, errors="ignore") as fh:
            for ln in fh:
                ln = ln.rstrip("\n")
                if not ln or ln.startswith("#"):
                    continue
                f = ln.split("\t")
                if f and f[0] == "replicon":          # cabecalho (pode repetir)
                    idx = {c: i for i, c in enumerate(f)}
                    continue
                if idx is None:
                    continue
                i = idx.get("hit_id")
                if i is not None and i < len(f) and f[i].strip():
                    hits.add(f[i].strip())
        if hits:
            break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="CONJScan/Chromosome CONJScan/Plasmids",
                    help="conjuntos de modelos (separados por espaco); CONJScan 2.x "
                         "tem subconjuntos Chromosome e Plasmids")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--log", default="")
    args = ap.parse_args()

    logfh = open(args.log, "a") if args.log else sys.stderr
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    def finish(hits):
        with open(args.out, "w") as out:
            out.write("contig\tstart\tend\n")
            for c, s, e in sorted(set(hits)):
                out.write(f"{c}\t{s}\t{e}\n")
        log(f"[conjscan_genome] hits escritos: {len(set(hits))} -> {args.out}", logfh)

    if not os.path.exists(args.assembly) or os.path.getsize(args.assembly) == 0:
        log("[conjscan_genome] sem assembly -> vazio. NAO e erro.", logfh)
        finish([]); return 0
    if not shutil.which("macsyfinder") or not shutil.which("prodigal"):
        log("[conjscan_genome] macsyfinder/prodigal ausentes -> vazio. NAO e erro.", logfh)
        finish([]); return 0

    work = args.workdir or os.path.join(os.path.dirname(args.out), "conjscan_genome_work")
    if os.path.isdir(work):
        shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    faa = os.path.join(work, "genome.faa")
    gff = os.path.join(work, "genome.gff")

    if not run_prodigal_all(args.assembly, faa, gff, logfh):
        log("[conjscan_genome] prodigal nao gerou proteinas -> vazio.", logfh)
        finish([]); return 0
    coords = faa_coords(faa)

    # 1 macsyfinder por conjunto de modelos (Chromosome + Plasmids); UNIAO dos
    # genes acertados. db-type unordered = recomendado p/ assembly multi-contig.
    all_hits = set()
    for ms in args.models.split():
        out_dir = os.path.join(work, "macsy_" + ms.replace("/", "_"))
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir, ignore_errors=True)
        cmd = ["macsyfinder", "--models", ms, "all",
               "--db-type", "unordered", "--sequence-db", faa,
               "-o", out_dir, "-w", "1"]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            log(f"[conjscan_genome] {ms}: macsyfinder retornou erro ({e}); seguindo.", logfh)
            continue
        h = parse_macsy_hits(out_dir)
        log(f"[conjscan_genome] {ms}: {len(h)} hit-genes", logfh)
        all_hits |= h

    hits = [coords[h] for h in all_hits if h in coords]
    log(f"[conjscan_genome] hit-genes de sistemas conjugativos: brutos={len(all_hits)} "
        f"mapeados={len(hits)}", logfh)
    finish(hits)
    if not args.workdir:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
