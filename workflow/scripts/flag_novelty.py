#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Flag de CANDIDATO A ESPÉCIE NOVA por vOTU — honesto e conservador.

ENCODIFICA A TRAVA DE COMPLETENESS: a definição operacional de espécie viral
(95% ANI sobre >=85% AF, padrão MIUViG) só é CONFIÁVEL quando o genoma está
completo o bastante. Em fragmento (low/medium), "não bater em referência" pode
ser apenas falta da parte conservada do genoma — então o status vira
'indeterminate_fragment', NUNCA 'candidate_novel_species'.

Critérios (todos derivados de saídas que o pipeline já produz):
  quality_ok      = CheckV miuvig_quality em {High-quality, Complete}
                    (ou completeness >= --min-completeness)
  is_phage        = portão de fago (mesma taxonomia geNomad/VITAP do enriquecimento)
  has_known_taxon = taxonomia (VITAP) alcança gênero/espécie (>= --known-rank ranks)

Status:
  not quality_ok                          -> indeterminate_fragment
  has_known_taxon                         -> classified (tem táxon de gênero/espécie)
  quality_ok & is_phage & not known_taxon -> candidate_novel_species
  quality_ok & not is_phage               -> non_phage_virus
Também reporta n_members (tamanho do cluster vOTU) como sinal de PREVALÊNCIA —
'recurrent' se n_members > 1 (aparece em >1 unidade/amostra).
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from build_gene_table import (classify_phage, DEFAULT_PHAGE_CLASSES,
                                   DEFAULT_EUK_CLASSES)
    import re as _re
except Exception:  # fallback defensivo
    DEFAULT_PHAGE_CLASSES = "Caudoviricetes"
    DEFAULT_EUK_CLASSES = "Herviviricetes|Megaviricetes|Revtraviricetes"
    import re as _re

    def classify_phage(tax, p, e, lenient):
        if not tax:
            return True
        if p.search(tax):
            return True
        if e.search(tax):
            return False
        return lenient

QUAL_OK = {"high-quality", "complete"}


def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, errors="ignore", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def taxon_depth(lineage):
    """Nº de ranks NÃO-vazios na taxonomia (ex.: geNomad
    'Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;' = 5).
    Também conta marcadores 'genus:'/'species:' do PhaGCN/VITAP."""
    if not lineage or lineage.strip().lower() in ("", "na", "unclassified", "-"):
        return 0
    parts = [p for p in _re.split(r"[;,]", lineage) if p.strip()
             and p.strip().lower() not in ("viruses", "na")]
    return len(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-bin", required=True, help="viral_bins_table.tsv")
    ap.add_argument("--checkv-bins", nargs="*", default=[], help="quality_summary.tsv do CheckV (bins)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--known-rank", type=int, default=6,
                    help="nº de ranks p/ considerar 'classificado a gênero+' (realm..genus=6)")
    ap.add_argument("--phage-classes", default=DEFAULT_PHAGE_CLASSES)
    ap.add_argument("--eukaryotic-classes", default=DEFAULT_EUK_CLASSES)
    args = ap.parse_args()

    phage_re = _re.compile(args.phage_classes, _re.I)
    euk_re = _re.compile(args.eukaryotic_classes, _re.I)

    # CheckV bins -> {contig_id: (miuvig_quality, completeness, completeness_method)}
    cv = {}
    for f in args.checkv_bins:
        for r in read_tsv(f):
            cid = r.get("contig_id", "")
            if cid:
                cv[cid] = (r.get("miuvig_quality", ""), r.get("completeness", ""),
                           r.get("completeness_method", ""))

    rows = read_tsv(args.per_bin)
    cols = ["votu_rep", "rep_unit", "checkv_quality", "miuvig_quality", "completeness",
            "is_phage", "taxonomy", "taxon_depth", "has_genus_level", "n_members",
            "prevalence", "novelty_status", "rationale"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    counts = {}
    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for b in rows:
            rep = b.get("votu_rep", "")
            unit = b.get("rep_unit", "")
            ckq = b.get("checkv_quality", "")
            miuvig, compl_cv, _meth = cv.get(unit, ("", "", ""))
            compl = b.get("completeness", "") or compl_cv
            try:
                compl_f = float(compl)
            except (TypeError, ValueError):
                compl_f = -1.0
            tax = b.get("taxonomy", "") or ""
            depth = taxon_depth(tax)
            try:
                n_members = int(b.get("n_members", 1) or 1)
            except ValueError:
                n_members = 1

            quality_ok = (miuvig.strip().lower() in QUAL_OK
                          or ckq.strip().lower() in QUAL_OK
                          or compl_f >= args.min_completeness)
            is_phage = classify_phage(tax, phage_re, euk_re, True)
            has_genus = depth >= args.known_rank

            if not quality_ok:
                status = "indeterminate_fragment"
                why = (f"completeness baixa (miuvig={miuvig or ckq or 'NA'}, "
                       f"compl={compl or 'NA'}); 95%ANI/85%AF não confiável em fragmento")
            elif has_genus:
                status = "classified"
                why = f"taxonomia alcança gênero/espécie (depth={depth})"
            elif not is_phage:
                status = "non_phage_virus"
                why = "qualidade OK mas taxonomia indica vírus eucariótico/não-fago"
            else:
                status = "candidate_novel_species"
                why = (f"High-quality+ ({miuvig or ckq}), fago, sem táxon de gênero/espécie "
                       f"(depth={depth}) -> novo em relação às referências")

            prevalence = "recurrent" if n_members > 1 else "singleton"
            counts[status] = counts.get(status, 0) + 1
            w.writerow({"votu_rep": rep, "rep_unit": unit, "checkv_quality": ckq,
                        "miuvig_quality": miuvig, "completeness": compl,
                        "is_phage": "Yes" if is_phage else "No", "taxonomy": tax or "-",
                        "taxon_depth": depth, "has_genus_level": "Yes" if has_genus else "No",
                        "n_members": n_members, "prevalence": prevalence,
                        "novelty_status": status, "rationale": why})

    summ = " | ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"[flag_novelty] {len(rows)} vOTUs -> {args.out}  [{summ}]")


if __name__ == "__main__":
    main()
