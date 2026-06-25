#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Converte a saida do IntegronFinder (rodado no GENOMA INTEIRO) em um TSV de
HITS por coordenada (contig, start, end) — usado por build_gene_table.py
(--integron-hits) para rotular genes como is_integron no teste de enriquecimento.

Cada integron vira 1 hit = span (min pos_beg, max pos_end) dos seus elementos.
Defensivo: dir ausente/sem .integrons -> TSV so com header.
"""
import argparse
import os


def parse_dir(integron_dir):
    hits = []
    if not integron_dir or not os.path.isdir(integron_dir):
        return hits
    for root, _d, files in os.walk(integron_dir):
        for fn in files:
            if not fn.endswith(".integrons"):
                continue
            fp = os.path.join(root, fn)
            try:
                with open(fp, errors="ignore") as fh:
                    header, idx, spans = None, {}, {}
                    for ln in fh:
                        ln = ln.rstrip("\n")
                        if not ln or ln.startswith("#"):
                            continue
                        f = ln.split("\t")
                        if header is None and ("ID_replicon" in f or "pos_beg" in f):
                            header = f
                            idx = {c: i for i, c in enumerate(f)}
                            continue
                        if header is None:
                            continue

                        def gv(col):
                            i = idx.get(col)
                            return f[i] if i is not None and i < len(f) else ""
                        rep = gv("ID_replicon")
                        iid = gv("ID_integron") or "1"
                        try:
                            pb, pe = int(float(gv("pos_beg"))), int(float(gv("pos_end")))
                        except ValueError:
                            continue
                        spans.setdefault((rep, iid), []).extend([pb, pe])
                    for (rep, _iid), pts in spans.items():
                        if rep and pts:
                            hits.append((rep, min(pts), max(pts)))
            except OSError:
                continue
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--integron-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    hits = parse_dir(args.integron_dir)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("contig\tstart\tend\n")
        for c, s, e in hits:
            out.write(f"{c}\t{s}\t{e}\n")
    print(f"[integron_to_hits] {len(hits)} integrons -> {args.out}")


if __name__ == "__main__":
    main()
