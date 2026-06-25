#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Filtra os contigs virais do geNomad por virus_score e reescreve o FASTA."""
import argparse


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fna", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--min-score", type=float, default=0.7)
    ap.add_argument("--out-fna", required=True)
    ap.add_argument("--out-ids", required=True)
    args = ap.parse_args()

    # virus_summary.tsv: col 'seq_name' e 'virus_score'
    keep = set()
    with open(args.summary) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = header.index("seq_name")
            i_score = header.index("virus_score")
        except ValueError:
            # se o cabeçalho mudar, mantém tudo
            i_name, i_score = 0, None
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            if i_score is None:
                keep.add(parts[i_name])
            else:
                try:
                    if float(parts[i_score]) >= args.min_score:
                        keep.add(parts[i_name])
                except (ValueError, IndexError):
                    keep.add(parts[i_name])

    n = 0
    with open(args.out_fna, "w") as out, open(args.out_ids, "w") as ids:
        for name, seq in read_fasta(args.fna):
            if not keep or name in keep:
                out.write(f">{name}\n{seq}\n")
                ids.write(name + "\n")
                n += 1
    print(f"[filter_genomad] mantidos {n} contigs virais (min_score={args.min_score})")


if __name__ == "__main__":
    main()
