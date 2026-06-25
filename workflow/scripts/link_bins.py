#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Une (link) os scaffolds de cada bin numa UNICA sequencia, separados por um
espacador de N's, para que o CheckV avalie o BIN INTEIRO (completude/qualidade
do vMAG) em vez de cada scaffold isolado.

Metodologia recomendada no artigo do vRhyme (Kieft et al. 2022, NAR):
  - junta os scaffolds do bin com um espacador de N's (default 1000)
  - o run de N's quebra ORFs na juncao (prodigal-gv, usado internamente pelo
    CheckV, nao chama genes atravessando N's), evitando genes quimericos
  - assim o CheckV estima a completude considerando todos os scaffolds juntos,
    ignorando a fragmentacao do bin.

Entrada : pasta com 1 .fna por bin (saida do collect_bins) + listing dos bins.
Saida   : 1 registro FASTA por bin  (>bin_id), pronto para o CheckV.

Bins de 1 unico scaffold saem sem espacador (identico ao scaffold).
"""
import argparse
import glob
import os


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


def bins_from_listing(listing):
    out = []
    if not os.path.exists(listing):
        return out
    with open(listing) as fh:
        next(fh, None)  # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                out.append(p[0])  # bin_id
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bindir", required=True,
                    help="pasta com <bin_id>.fna (collect_bins)")
    ap.add_argument("--listing", required=True,
                    help="TSV dos bins (collect_bins) p/ saber a ordem/ids")
    ap.add_argument("--out-fna", required=True)
    ap.add_argument("--spacer-n", type=int, default=1000,
                    help="num. de N's entre scaffolds (default 1000)")
    args = ap.parse_args()

    spacer = "N" * args.spacer_n
    os.makedirs(os.path.dirname(args.out_fna) or ".", exist_ok=True)

    bin_ids = bins_from_listing(args.listing)
    # fallback: se nao houver listing, varre a pasta
    if not bin_ids:
        bin_ids = [os.path.basename(f)[:-4]
                   for f in sorted(glob.glob(os.path.join(args.bindir, "*.fna")))]

    n_bins, n_multi = 0, 0
    with open(args.out_fna, "w") as out:
        for bin_id in bin_ids:
            f = os.path.join(args.bindir, bin_id + ".fna")
            scaffolds = [s for _, s in read_fasta(f)]
            if not scaffolds:
                continue
            if len(scaffolds) > 1:
                n_multi += 1
            linked = spacer.join(scaffolds)
            out.write(f">{bin_id}\n{linked}\n")
            n_bins += 1

    print(f"[link_bins] {n_bins} bins linkados "
          f"({n_multi} multi-scaffold, espacador={args.spacer_n} N's) -> {args.out_fna}")


if __name__ == "__main__":
    main()
