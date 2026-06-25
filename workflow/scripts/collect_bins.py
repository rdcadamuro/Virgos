#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Coleta bins de vRhyme e VAMB numa pasta única por amostra, renomeando
contigs e arquivos para rastrear origem: <sample>__<binner>__bin<N>.

Gera um TSV: bin_id, sample, binner, source_file, n_contigs.
Cada contig é renomeado para <bin_id>__<contig_original> para que o CheckV
e os passos seguintes saibam a qual bin cada sequência pertence.
"""
import argparse
import glob
import os


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


def collect(binner, indir, sample, outdir, rows):
    if not indir or not os.path.isdir(indir):
        return
    files = sorted(glob.glob(os.path.join(indir, "*.fna")) +
                   glob.glob(os.path.join(indir, "*.fasta")) +
                   glob.glob(os.path.join(indir, "*.fa")))
    for k, f in enumerate(files, 1):
        bin_id = f"{sample}__{binner}__bin{k}"
        out_path = os.path.join(outdir, bin_id + ".fna")
        n = 0
        with open(out_path, "w") as out:
            for name, seq in read_fasta(f):
                out.write(f">{bin_id}__{name}\n{seq}\n")
                n += 1
        if n == 0:
            os.remove(out_path)
            continue
        rows.append((bin_id, sample, binner, os.path.basename(f), n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vrhyme-dir", required=True)
    ap.add_argument("--vamb-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--listing", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rows = []
    collect("vrhyme", args.vrhyme_dir, args.sample, args.out_dir, rows)
    collect("vamb", args.vamb_dir, args.sample, args.out_dir, rows)

    with open(args.listing, "w") as out:
        out.write("bin_id\tsample\tbinner\tsource_file\tn_contigs\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    print(f"[collect_bins] {args.sample}: {len(rows)} bins coletados")


if __name__ == "__main__":
    main()
