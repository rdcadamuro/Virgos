#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tabela de GENES do genoma completo, rotulada para o teste de enriquecimento.

Universo = TODOS os genes do genoma (geNomad annotate/_genes.tsv).
Rótulos por gene:
  - is_viral : gene cai numa região viral (contig totalmente viral OU dentro de
               um provírus) — geNomad virus_summary + provirus.tsv
  - is_arg   : geNomad annotation_amr preenchido  [+ opcional: hit ABRicate no contig]
  - is_mge   : geNomad annotation_conjscan preenchido OU descrição casa MGE
               [+ opcional: gene dentro de uma mobile_element_region]

geNomad já anota AMR/conjugação para o GENOMA INTEIRO (viral + não-viral) de
forma consistente — é o background "de graça". ABRicate/regiões de mobiloma no
genoma inteiro são augmentações opcionais (mesmo critério aplicado a todos os genes).

Saída: TSV (gene, sample, contig, start, end, is_viral, is_arg, is_mge) — entrada
do permutation_enrichment.py.
"""
import argparse
import os
import re

NA = {"", "NA", "na", "-", "nan", "None"}


def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


def gene_contig(gene_id):
    return gene_id.rsplit("_", 1)[0] if "_" in gene_id else gene_id


def parse_viral_regions(virus_summary, provirus):
    """Retorna (set de contigs totalmente virais, lista de (src, start, end) provírus)."""
    full, prov = set(), []
    for r in read_tsv(virus_summary):
        sn = r.get("seq_name", "")
        if "|provirus_" in sn:
            # ex.: NODE_33...|provirus_131_14917
            m = re.search(r"\|provirus_(\d+)_(\d+)$", sn)
            src = sn.split("|provirus_")[0]
            if m:
                prov.append((src, int(m.group(1)), int(m.group(2))))
        elif sn:
            full.add(sn)
    # provirus.tsv (source_seq, start, end) — reforça
    for r in read_tsv(provirus):
        src = r.get("source_seq") or r.get("seq_name", "")
        try:
            s, e = int(r.get("start", 0)), int(r.get("end", 0))
            if src and e > 0:
                prov.append((src.split("|provirus_")[0], s, e))
        except ValueError:
            pass
    return full, prov


def load_coord_hits(path, seq_key_candidates=("sequence", "seq", "contig")):
    """Retorna {contig: [(start,end), ...]} de uma tabela com coordenadas."""
    out = {}
    rows = read_tsv(path)
    if not rows:
        return out
    cols = rows[0].keys()
    sk = next((c for c in cols if c.lower() in seq_key_candidates), None)
    if sk is None:
        return out
    for r in rows:
        try:
            s, e = int(float(r.get("start", 0))), int(float(r.get("end", 0)))
        except ValueError:
            continue
        out.setdefault(r.get(sk, ""), []).append((min(s, e), max(s, e)))
    return out


def overlaps(gs, ge, intervals):
    return any(not (ge < s or gs > e) for s, e in intervals)


# Regex por CATEGORIA (sobre annotation_description do geNomad, case-insensitive)
_CAT_RE = {
    "is_integrase":  re.compile(r"integrase", re.I),
    "is_transposase": re.compile(r"transposase|insertion sequence|\bIS\d", re.I),
    "is_recombinase": re.compile(r"recombinase|relaxase|resolvase", re.I),
}

# Colunas de categoria emitidas (ordem fixa p/ a tabela geral ficar estavel).
# geNomad cobre o GENOMA INTEIRO p/: arg, conjugation, integrase, transposase, recombinase.
# integron / conjugative_system / IS dependem de ferramentas no genoma inteiro
# (--integron-hits / --conjscan-hits / --is-hits); sem elas ficam 0 (-> UNDERPOWERED).
_CATEGORIES = ["is_arg", "is_conjugation", "is_integrase", "is_transposase",
               "is_recombinase", "is_integron", "is_conjugative_system",
               "is_IS", "is_mge"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="geNomad annotate/<s>_genes.tsv")
    ap.add_argument("--virus-summary", required=True)
    ap.add_argument("--provirus", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mge-keywords",
                    default="integrase|transposase|recombinase|relaxase|resolvase|mobiliz|conjugat|insertion sequence")
    ap.add_argument("--abricate", default="", help="(opcional) ABRicate dedup no genoma (coords) -> is_arg")
    ap.add_argument("--mge-regions", default="", help="(opcional) mobile_element_regions no genoma -> is_mge")
    ap.add_argument("--integron-hits", default="", help="(opcional) IntegronFinder no genoma (coords) -> is_integron")
    ap.add_argument("--conjscan-hits", default="", help="(opcional) CONJScan no genoma (coords) -> is_conjugative_system")
    ap.add_argument("--is-hits", default="", help="(opcional) MobileElementFinder no genoma (coords) -> is_IS")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mge_re = re.compile(args.mge_keywords, re.IGNORECASE)
    full_viral, prov = parse_viral_regions(args.virus_summary, args.provirus)
    abr = load_coord_hits(args.abricate) if args.abricate else {}
    mge_reg = load_coord_hits(args.mge_regions) if args.mge_regions else {}
    int_hits = load_coord_hits(args.integron_hits) if args.integron_hits else {}
    conj_hits = load_coord_hits(args.conjscan_hits) if args.conjscan_hits else {}
    is_hits = load_coord_hits(args.is_hits) if args.is_hits else {}

    genes = read_tsv(args.genes)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    totals = {c: 0 for c in _CATEGORIES}
    n_viral = 0
    header = ["gene", "sample", "contig", "start", "end", "is_viral"] + _CATEGORIES
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for g in genes:
            gid = g.get("gene", "")
            contig = gene_contig(gid)
            try:
                gs, ge = int(g.get("start", 0)), int(g.get("end", 0))
            except ValueError:
                gs, ge = 0, 0
            desc = g.get("annotation_description", "") or ""

            viral = contig in full_viral or any(
                src == contig and not (ge < s or gs > e) for src, s, e in prov)

            # --- categorias do geNomad (genoma inteiro) ---
            arg = g.get("annotation_amr", "") not in NA
            if not arg and contig in abr:
                arg = overlaps(gs, ge, abr[contig])
            conj = g.get("annotation_conjscan", "") not in NA
            integrase = bool(_CAT_RE["is_integrase"].search(desc))
            transposase = bool(_CAT_RE["is_transposase"].search(desc))
            recombinase = bool(_CAT_RE["is_recombinase"].search(desc))

            # --- categorias de ferramentas no genoma inteiro (coords, opcional) ---
            integron = overlaps(gs, ge, int_hits.get(contig, []))
            conj_sys = overlaps(gs, ge, conj_hits.get(contig, []))
            is_elem = overlaps(gs, ge, is_hits.get(contig, []))

            # --- MGE (uniao) ---
            mge = (conj or integrase or transposase or recombinase
                   or integron or conj_sys or is_elem
                   or bool(mge_re.search(desc)))
            if not mge and contig in mge_reg:
                mge = overlaps(gs, ge, mge_reg[contig])

            row = {
                "is_arg": arg, "is_conjugation": conj, "is_integrase": integrase,
                "is_transposase": transposase, "is_recombinase": recombinase,
                "is_integron": integron, "is_conjugative_system": conj_sys,
                "is_IS": is_elem, "is_mge": mge,
            }
            n_viral += viral
            for c in _CATEGORIES:
                totals[c] += int(row[c])
            out.write(f"{gid}\t{args.sample}\t{contig}\t{gs}\t{ge}\t{int(viral)}\t"
                      + "\t".join(str(int(row[c])) for c in _CATEGORIES) + "\n")

    cat_str = " ".join(f"{c.replace('is_', '')}={totals[c]}" for c in _CATEGORIES)
    print(f"[build_gene_table] {args.sample}: {len(genes)} genes "
          f"(viral={n_viral}; {cat_str}) -> {args.out}")


if __name__ == "__main__":
    main()
