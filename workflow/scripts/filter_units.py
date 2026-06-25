#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define as UNIDADES virais finais aplicando o filtro do CheckV
(>= min_viral_genes genes virais / hallmark).

  --mode metagenome : unidade = BIN; soma viral_genes dos contigs do bin.
  --mode genome     : unidade = CONTIG viral com >= min_viral_genes.

Saídas:
  passed.fna       : sequências aprovadas
  passed_units.tsv : unit_id, sample, binner, viral_genes, completeness, checkv_quality
"""
import argparse
import os

QUALITY_RANK = {"Complete": 5, "High-quality": 4, "Medium-quality": 3,
                "Low-quality": 2, "Not-determined": 1, "": 0}


def read_fasta(path):
    name, seq = None, []
    if not os.path.exists(path):
        return
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


def load_checkv(path):
    """contig_id -> dict(viral_genes, completeness, quality)."""
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0]:
                continue
            cid = p[idx.get("contig_id", 0)]
            try:
                vg = int(float(p[idx["viral_genes"]])) if "viral_genes" in idx else 0
            except (ValueError, KeyError):
                vg = 0
            try:
                comp = float(p[idx["completeness"]]) if p[idx["completeness"]] not in ("", "NA") else 0.0
            except (ValueError, KeyError):
                comp = 0.0
            qual = p[idx["checkv_quality"]] if "checkv_quality" in idx else ""
            out[cid] = {"viral_genes": vg, "completeness": comp, "quality": qual}
    return out


def write_outputs(units, fnas, out_fna, out_units):
    with open(out_units, "w") as u:
        u.write("unit_id\tsample\tbinner\tviral_genes\tcompleteness\tcheckv_quality\n")
        for r in units:
            u.write("\t".join(map(str, r)) + "\n")
    with open(out_fna, "w") as out:
        for blob in fnas:
            out.write(blob)


def run_metagenome(args, checkv):
    """CheckV foi rodado no BIN INTEIRO (scaffolds unidos por N's, link_bins).
    Logo a tabela do CheckV tem contig_id == bin_id e a completude/qualidade
    ja sao do bin como unidade (vMAG) — sem agregacao por contig.
    As sequencias gravadas em passed.fna continuam sendo os scaffolds REAIS do
    bin (sem os N's), para derep / PhaBox2 / ViTAP / abundancia.
    """
    bins = {}  # bin_id -> binner
    with open(args.listing) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                bins[p[0]] = p[2]

    units, fnas = [], []
    for bin_id, binner in bins.items():
        f = os.path.join(args.bindir, bin_id + ".fna")
        if not os.path.exists(f):
            continue
        contigs = list(read_fasta(f))
        # CheckV do bin inteiro: chave = bin_id (vide link_bins.py)
        cv = checkv.get(bin_id, {"viral_genes": 0, "completeness": 0.0, "quality": ""})
        if cv["viral_genes"] >= args.min_viral_genes:
            units.append((bin_id, args_sample_from_bin(bin_id), binner,
                          cv["viral_genes"], f"{cv['completeness']:.1f}", cv["quality"]))
            fnas.append("".join(f">{n}\n{s}\n" for n, s in contigs))
    return units, fnas


def args_sample_from_bin(bin_id):
    # bin_id = <sample>__<binner>__bin<N>
    return bin_id.split("__")[0]


def run_genome(args, checkv):
    units, fnas = [], []
    for name, seq in read_fasta(args.viral_fna):
        cv = checkv.get(name, {"viral_genes": 0, "completeness": 0.0, "quality": ""})
        if cv["viral_genes"] >= args.min_viral_genes:
            units.append((name, args.sample, "genomad-contig",
                          cv["viral_genes"], f"{cv['completeness']:.1f}", cv["quality"]))
            fnas.append(f">{name}\n{seq}\n")
    return units, fnas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["metagenome", "genome"], required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--min-viral-genes", type=int, default=1)
    ap.add_argument("--out-fna", required=True)
    ap.add_argument("--out-units", required=True)
    # metagenome
    ap.add_argument("--listing")
    ap.add_argument("--bindir")
    # genome
    ap.add_argument("--viral-fna")
    ap.add_argument("--sample")
    args = ap.parse_args()

    checkv = load_checkv(args.summary)
    if args.mode == "metagenome":
        units, fnas = run_metagenome(args, checkv)
    else:
        units, fnas = run_genome(args, checkv)

    write_outputs(units, fnas, args.out_fna, args.out_units)
    print(f"[filter_units:{args.mode}] {len(units)} unidades aprovadas "
          f"(min_viral_genes={args.min_viral_genes})")


if __name__ == "__main__":
    main()
