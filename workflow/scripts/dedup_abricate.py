#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deduplica hits do ABRicate por LOCUS + CHECKPOINT auto-verificavel.

Ao rodar VARIOS bancos do ABRicate (NCBI, CARD, ResFinder, ARG-ANNOT, MEGARES,
VFDB, Victors, PlasmidFinder, BacMet2, ...), o MESMO gene fisico costuma ser
detectado por mais de um banco -> linhas redundantes. Este script:

  1) agrupa hits da MESMA sequencia por SOBREPOSICAO RECIPROCA (cada intervalo
     precisa cobrir o outro em >= --overlap) via UNION-FIND. Isso:
       - funde "mesmo gene, bancos diferentes" (coordenadas ~iguais);
       - NAO encadeia genes adjacentes distintos (A~B, B~C mas A!~C) nem funde
         um gene pequeno aninhado num hit maior (cobertura reciproca baixa).
     -> corrige o encadeamento transitivo do clustering guloso antigo.
  2) mantem UMA linha por locus: o melhor hit (maior %id x %cov), listando todos
     os bancos e nomes de gene que concordaram (concordancia = confianca).
  3) grava um CHECKPOINT (JSON) com contagens (bruto -> loci, por banco) e uma
     VERIFICACAO: nenhuma sobreposicao residual na saida (status=PASS sse 0).
  4) (opcional) grava um RESUMO POR BANCO (--out-by-db): quantos hits cada banco
     trouxe e em quantos loci nao-redundantes ele participa.

Entrada : 1+ tabelas do ABRicate concatenadas (com headers '#FILE...' repetidos).
Saida   : tabela deduplicada por locus (+ checkpoint + resumo por banco).
"""
import argparse
import csv
import json
import os
from collections import defaultdict

# Colunas padrao do ABRicate (fallback se nao houver linha de header)
COLS = ["FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
        "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE",
        "ACCESSION", "PRODUCT", "RESISTANCE"]


def read_hits(paths):
    """Le 1+ tabelas do ABRicate. Trata headers '#FILE...' repetidos (1 por banco)."""
    hits = []
    for p in paths:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        with open(p, newline="", errors="ignore") as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    header = [c.lstrip("#") for c in line.lstrip("#").split("\t")]
                    continue
                parts = line.split("\t")
                cols = header if header else COLS
                row = dict(zip(cols, parts))
                try:
                    row["_start"] = int(row.get("START", 0))
                    row["_end"] = int(row.get("END", 0))
                    row["_id"] = float(row.get("%IDENTITY", 0) or 0)
                    row["_cov"] = float(row.get("%COVERAGE", 0) or 0)
                except ValueError:
                    continue
                if row["_end"] < row["_start"]:
                    row["_start"], row["_end"] = row["_end"], row["_start"]
                hits.append(row)
    return hits


def reciprocal_overlap(a, b):
    """min das duas fracoes de cobertura (0..1). Alto so se os spans sao parecidos."""
    s = max(a["_start"], b["_start"])
    e = min(a["_end"], b["_end"])
    ov = max(0, e - s + 1)
    if ov <= 0:
        return 0.0
    la = max(1, a["_end"] - a["_start"] + 1)
    lb = max(1, b["_end"] - b["_start"] + 1)
    return min(ov / la, ov / lb)


def cluster_locus(hits, overlap):
    """Agrupa hits da MESMA sequencia por sobreposicao reciproca (union-find).
    Independe de ordem; nao encadeia (A~B,B~C mas A!~C ficam separados se A!~C
    e o union-find so une pares que de fato se sobrepoem reciprocamente)."""
    by_seq = defaultdict(list)
    for h in hits:
        by_seq[h.get("SEQUENCE", "")].append(h)
    clusters = []
    for _seq, hs in by_seq.items():
        n = len(hs)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            parent[find(x)] = find(y)

        for i in range(n):
            for j in range(i + 1, n):
                if reciprocal_overlap(hs[i], hs[j]) >= overlap:
                    union(i, j)
        groups = defaultdict(list)
        for i in range(n):
            groups[find(i)].append(hs[i])
        clusters.extend(groups.values())
    return clusters


def verify_no_residual(loci, overlap):
    """Auto-verificacao: nenhum par de loci da MESMA sequencia pode se sobrepor
    reciprocamente >= overlap (se sobrepoe, a clusterizacao falhou)."""
    by_seq = defaultdict(list)
    for lc in loci:
        by_seq[lc["sequence"]].append(lc)
    residual = 0
    for _seq, ls in by_seq.items():
        for i in range(len(ls)):
            for j in range(i + 1, len(ls)):
                a = {"_start": ls[i]["start"], "_end": ls[i]["end"]}
                b = {"_start": ls[j]["start"], "_end": ls[j]["end"]}
                if reciprocal_overlap(a, b) >= overlap:
                    residual += 1
    return residual


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inputs", nargs="+", required=True)
    ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-checkpoint", default="", help="JSON com contagens + verificacao")
    ap.add_argument("--out-by-db", default="", help="TSV: contribuicao por banco")
    ap.add_argument("--dbs-run", default="", help="lista (espaco) de bancos executados; "
                    "incluidos no resumo mesmo com 0 hits (proveniencia: o que foi comparado)")
    args = ap.parse_args()

    hits = read_hits(args.inputs)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # contagem bruta por banco
    raw_by_db = defaultdict(int)
    for h in hits:
        raw_by_db[h.get("DATABASE", "?")] += 1

    out_cols = ["sequence", "start", "end", "strand", "best_gene", "pct_identity",
                "pct_coverage", "n_databases", "databases", "all_genes",
                "product", "resistance"]
    loci = []
    loci_by_db = defaultdict(int)  # em quantos loci cada banco participa
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(out_cols)
        for cl in cluster_locus(hits, args.overlap):
            best = max(cl, key=lambda x: x["_id"] * x["_cov"])
            dbs = sorted({h.get("DATABASE", "") for h in cl if h.get("DATABASE")})
            genes = sorted({h.get("GENE", "") for h in cl if h.get("GENE")})
            start = min(h["_start"] for h in cl)
            end = max(h["_end"] for h in cl)
            for db in dbs:
                loci_by_db[db] += 1
            w.writerow([best.get("SEQUENCE", ""), start, end, best.get("STRAND", ""),
                        best.get("GENE", ""), f"{best['_id']:.1f}", f"{best['_cov']:.1f}",
                        len(dbs), ";".join(dbs), ";".join(genes),
                        best.get("PRODUCT", ""), best.get("RESISTANCE", "")])
            loci.append({"sequence": best.get("SEQUENCE", ""), "start": start, "end": end})

    residual = verify_no_residual(loci, args.overlap)
    status = "PASS" if residual == 0 else "FAIL"

    # bancos executados (proveniencia): inclui os de --dbs-run mesmo com 0 hits
    dbs_run = args.dbs_run.split() if args.dbs_run else []
    all_dbs = sorted(set(raw_by_db) | set(loci_by_db) | set(dbs_run))

    if args.out_by_db:
        os.makedirs(os.path.dirname(args.out_by_db) or ".", exist_ok=True)
        with open(args.out_by_db, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["database", "n_raw_hits", "n_loci_supported"])
            for db in all_dbs:
                w.writerow([db, raw_by_db.get(db, 0), loci_by_db.get(db, 0)])

    if args.out_checkpoint:
        os.makedirs(os.path.dirname(args.out_checkpoint) or ".", exist_ok=True)
        cp = {
            "overlap_threshold": args.overlap,
            "n_databases_run": len(dbs_run) if dbs_run else len(all_dbs),
            "databases_run": dbs_run or all_dbs,
            "n_raw_hits": len(hits),
            "n_loci": len(loci),
            "raw_by_db": dict(raw_by_db),
            "loci_by_db": dict(loci_by_db),
            "verification": {"residual_overlaps": residual, "status": status},
        }
        with open(args.out_checkpoint, "w") as fh:
            json.dump(cp, fh, indent=2)

    print(f"[dedup_abricate] {len(hits)} hits -> {len(loci)} loci unicos "
          f"(overlap>={args.overlap}) | verificacao={status} (residual={residual})")


if __name__ == "__main__":
    main()
