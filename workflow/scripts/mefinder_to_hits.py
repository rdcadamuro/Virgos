#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MobileElementFinder .csv (rodado no GENOMA INTEIRO) -> TSV de HITS por
coordenada (contig, start, end), usado por build_gene_table.py (--is-hits) para
rotular genes como is_IS no teste de enriquecimento.

Sem isto, a categoria is_IS fica estruturalmente vazia (0 -> UNDERPOWERED),
pois o mefinder so rodava nos vOTUs (sem background nao-viral, sem enriquecimento).

Defensivo: arquivo ausente/vazio/sem colunas de coordenada -> TSV so com header.
"""
import argparse
import csv
import os


def parse_csv(path):
    """mefinder.csv (linhas '#' = comentario do cabecalho MGEdb). Detecta as
    colunas de contig/start/end de forma tolerante (nomes variam por versao)."""
    hits = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return hits
    with open(path, errors="ignore") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return hits
    rows = list(csv.reader(lines))
    if not rows:
        return hits
    low = [h.strip().lower() for h in rows[0]]

    def col(*names):
        for nm in names:
            if nm in low:
                return low.index(nm)
        return None

    ci = col("contig")
    if ci is None:
        ci = next((i for i, h in enumerate(low) if "contig" in h or "sequence" in h), None)
    si = col("start", "contig_start", "begin")
    ei = col("end", "contig_end", "stop")
    if ci is None or si is None or ei is None:
        return hits
    for r in rows[1:]:
        if len(r) <= max(ci, si, ei) or not r[ci].strip():
            continue
        try:
            s, e = int(float(r[si])), int(float(r[ei]))
        except ValueError:
            continue
        hits.append((r[ci].strip(), min(s, e), max(s, e)))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    hits = parse_csv(args.csv)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("contig\tstart\tend\n")
        for c, s, e in hits:
            out.write(f"{c}\t{s}\t{e}\n")
    print(f"[mefinder_to_hits] {len(hits)} IS/transposons -> {args.out}")


if __name__ == "__main__":
    main()
