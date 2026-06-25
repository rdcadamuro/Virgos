#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Limpa a saída do CoverM: renomeia colunas '<path>/SAMPLE.votu.sorted.bam Method'
para apenas 'SAMPLE' e renomeia a 1ª coluna para 'votu'."""
import argparse
import os
import re


def clean_col(col):
    # remove o nome do método ao final (TPM, Mean, Covered Fraction ...)
    base = col
    base = re.sub(r"\s+(TPM|Mean|Covered Fraction|Relative Abundance \(%\)|Read Count|RPKM|Variance|Length)$",
                  "", base, flags=re.IGNORECASE)
    base = os.path.basename(base.strip())
    base = re.sub(r"\.votu\.sorted\.bam$", "", base)
    base = re.sub(r"\.sorted\.bam$", "", base)
    base = re.sub(r"\.bam$", "", base)
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if not os.path.exists(args.raw) or os.path.getsize(args.raw) == 0:
        with open(args.out, "w") as o:
            o.write("votu\n")
        return

    with open(args.raw) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        new_header = ["votu"] + [clean_col(c) for c in header[1:]]
        with open(args.out, "w") as o:
            o.write("\t".join(new_header) + "\n")
            for line in fh:
                o.write(line if line.endswith("\n") else line + "\n")
    print(f"[clean_coverm] {args.out}")


if __name__ == "__main__":
    main()
