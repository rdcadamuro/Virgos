#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dereplicação de sequências virais em vOTUs.

Implementa o procedimento padrão do CheckV (anicalc + aniclust):
  1) makeblastdb + blastn all-vs-all (megablast)
  2) calcula ANI e cobertura (alignment fraction) entre pares
  3) clusteriza por greedy centroid: ANI >= min_ani e AF >= min_af
O representante de cada cluster é a sequência mais longa.

Requer: blast (makeblastdb, blastn) no PATH.
Refs: Nayfach et al. 2021 (CheckV); scripts anicalc.py/aniclust.py.
"""
import argparse
import os
import subprocess
from collections import defaultdict


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def run(cmd):
    print("  $", " ".join(cmd))
    subprocess.run(cmd, check=True)


def compute_ani(blast_tsv, lengths):
    """Recria anicalc: para cada par (qseqid, sseqid), ANI ponderada por
    comprimento alinhado e AF = bases alinhadas / menor sequência."""
    # acumula por par
    pair = defaultdict(lambda: {"len": 0, "ident": 0.0})
    with open(blast_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            if q == s:
                continue
            pid = float(f[2])
            alen = int(f[3])
            key = (q, s)
            pair[key]["len"] += alen
            pair[key]["ident"] += pid * alen
    edges = []
    for (q, s), v in pair.items():
        if v["len"] == 0:
            continue
        ani = v["ident"] / v["len"]
        qlen = lengths.get(q, 1)
        slen = lengths.get(s, 1)
        af = 100.0 * v["len"] / min(qlen, slen)
        edges.append((q, s, ani, min(af, 100.0)))
    return edges


def cluster(seqs, edges, min_ani, min_af):
    """aniclust: greedy por centroides ordenados por comprimento (desc)."""
    lengths = {n: len(s) for n, s in seqs}
    order = sorted(lengths, key=lambda n: lengths[n], reverse=True)

    # adjacência de quem satisfaz os limiares
    nbr = defaultdict(set)
    for q, s, ani, af in edges:
        if ani >= min_ani and af >= min_af:
            nbr[q].add(s)
            nbr[s].add(q)

    assigned = set()
    clusters = []  # (rep, [members])
    for rep in order:
        if rep in assigned:
            continue
        members = [rep]
        assigned.add(rep)
        for m in nbr[rep]:
            if m not in assigned:
                members.append(m)
                assigned.add(m)
        clusters.append((rep, members))
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fna", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=85.0)
    ap.add_argument("--out-votus", required=True)
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-map", required=True)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    seqs = list(read_fasta(args.fna))
    lengths = {n: len(s) for n, s in seqs}

    db = os.path.join(args.workdir, "db")
    blast_out = os.path.join(args.workdir, "blast.tsv")

    run(["makeblastdb", "-in", args.fna, "-dbtype", "nucl", "-out", db])
    run(["blastn", "-query", args.fna, "-db", db,
         "-outfmt", "6 std", "-max_target_seqs", "10000",
         "-num_threads", str(args.threads),
         "-perc_identity", "90", "-out", blast_out])

    edges = compute_ani(blast_out, lengths)
    clusters = cluster(seqs, edges, args.min_ani, args.min_af)

    seqdict = dict(seqs)
    with open(args.out_votus, "w") as fv, \
         open(args.out_clusters, "w") as fc, \
         open(args.out_map, "w") as fm:
        fc.write("representative\tmembers\n")
        fm.write("contig\tvotu_rep\n")
        for rep, members in clusters:
            fv.write(f">{rep}\n{seqdict[rep]}\n")
            fc.write(f"{rep}\t{','.join(members)}\n")
            for m in members:
                fm.write(f"{m}\t{rep}\n")

    print(f"[dereplicate] {len(seqs)} seqs -> {len(clusters)} vOTUs "
          f"(ANI>={args.min_ani}, AF>={args.min_af})")


if __name__ == "__main__":
    main()
