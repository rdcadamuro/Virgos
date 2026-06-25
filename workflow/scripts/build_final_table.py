#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bloco final: tabela por vOTU + resumo por amostra.

Junta:
  - clusters.tsv (derep)            -> membros de cada vOTU
  - all_passed_units.tsv            -> sample, binner, qualidade, completeness por unidade
  - phatyp_prediction.tsv (PhaBox2) -> lifestyle temperate(lysogenic)/virulent(lytic)
  - vitap_taxonomy.tsv (ViTAP)      -> taxonomia
  - [metagenoma] votu_abundance_tpm.tsv -> prevalência e riqueza por amostra
"""
import argparse
import csv
import os
from collections import defaultdict


def sniff_read(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path, newline="") as fh:
        sample = fh.read(2048)
        fh.seek(0)
        delim = "\t" if "\t" in sample else ","
        rows = [r for r in csv.reader(fh, delimiter=delim) if r]
    if not rows:
        return [], []
    return rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]


def load_units(path):
    _, rows = sniff_read(path)
    out = {}
    for r in rows:
        out[r["unit_id"]] = {
            "sample": r.get("sample", ""),
            "binner": r.get("binner", ""),
            "checkv_quality": r.get("checkv_quality", ""),
            "completeness": r.get("completeness", ""),
        }
    return out


def unit_of(contig, unit_ids_sorted):
    return next((u for u in unit_ids_sorted
                 if contig == u or contig.startswith(u + "__")), None)


def load_phatyp(path):
    header, rows = sniff_read(path)
    out = {}
    if not header:
        return out
    id_col = next((c for c in header if c.lower() in
                   ("accession", "contig", "seqid", "id", "sequence")), header[0])
    type_col = next((c for c in header if c.lower() in
                     ("type", "pred", "prediction", "lifestyle", "phatyp")), None)
    score_col = next((c for c in header if "score" in c.lower()), None)
    for r in rows:
        raw = (r.get(type_col, "") if type_col else "").strip().lower()
        if raw.startswith("temp"):
            ls = "lysogenic"
        elif raw.startswith(("viru", "lyt")):
            ls = "lytic"
        else:
            ls = raw or "unknown"
        out[r.get(id_col, "")] = (ls, r.get(score_col, "") if score_col else "")
    return out


def load_vitap(path):
    header, rows = sniff_read(path)
    if not header:
        return {}
    id_col = next((c for c in header if c.lower() in
                   ("contig", "accession", "seqid", "id", "query")), header[0])
    tax_col = next((c for c in header if any(k in c.lower()
                    for k in ("tax", "lineage", "classif"))), None)
    out = {}
    for r in rows:
        cid = r.get(id_col, "")
        out[cid] = r.get(tax_col, "") if tax_col else \
            ";".join(v for k, v in r.items() if k != id_col)
    return out


def load_abundance(path):
    """votu -> {sample: tpm}; e lista de samples (colunas)."""
    header, rows = sniff_read(path)
    if not header or len(header) < 2:
        return {}, []
    samples = header[1:]
    out = {}
    for r in rows:
        votu = r.get(header[0], "")
        out[votu] = {s: float(r.get(s, 0) or 0) for s in samples}
    return out, samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--votu-map")
    ap.add_argument("--units", required=True)
    ap.add_argument("--phatyp", required=True)
    ap.add_argument("--vitap", required=True)
    ap.add_argument("--abundance-tpm", default=None)
    ap.add_argument("--out-per-bin", required=True)
    ap.add_argument("--out-per-sample", required=True)
    args = ap.parse_args()

    units = load_units(args.units)
    unit_ids_sorted = sorted(units, key=len, reverse=True)
    phatyp = load_phatyp(args.phatyp)
    vitap = load_vitap(args.vitap)
    abund, abund_samples = load_abundance(args.abundance_tpm) if args.abundance_tpm else ({}, [])

    _, crows = sniff_read(args.clusters)
    os.makedirs(os.path.dirname(args.out_per_bin), exist_ok=True)

    cols = ["votu_rep", "rep_unit", "rep_sample", "all_samples", "binner",
            "checkv_quality", "completeness", "lifestyle", "lifestyle_score",
            "taxonomy", "n_members"]
    if abund:
        cols += ["prevalence", "mean_tpm"]

    sample_life = defaultdict(lambda: defaultdict(int))
    sample_qual = defaultdict(lambda: defaultdict(int))
    sample_detected = defaultdict(int)   # riqueza: nº de vOTUs com TPM>0

    with open(args.out_per_bin, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for r in crows:
            rep = r.get("representative", "")
            members = [m for m in r.get("members", "").split(",") if m]
            rep_unit = unit_of(rep, unit_ids_sorted)
            info = units.get(rep_unit, {})
            rep_sample = info.get("sample", "")

            samples_here = sorted({
                units.get(unit_of(m, unit_ids_sorted), {}).get("sample", "")
                for m in members
            } - {""})

            ls, score = phatyp.get(rep, ("unknown", ""))
            tax = vitap.get(rep, "")
            row = [rep, rep_unit or "", rep_sample, ";".join(samples_here),
                   info.get("binner", ""), info.get("checkv_quality", ""),
                   info.get("completeness", ""), ls, score, tax, len(members)]

            if abund:
                tpms = abund.get(rep, {})
                present = [s for s, v in tpms.items() if v > 0]
                for s in present:
                    sample_detected[s] += 1
                mean_tpm = (sum(tpms.values()) / len(tpms)) if tpms else 0.0
                row += [len(present), f"{mean_tpm:.4f}"]

            w.writerow(row)

            if rep_sample:
                sample_life[rep_sample][ls] += 1
                sample_qual[rep_sample][info.get("checkv_quality", "NA") or "NA"] += 1

    qtiers = ["Complete", "High-quality", "Medium-quality",
              "Low-quality", "Not-determined", "NA"]
    with open(args.out_per_sample, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        head = ["sample", "n_votus_origin", "n_lysogenic", "n_lytic", "n_unknown"] \
            + [f"q_{t}" for t in qtiers]
        if abund:
            head.append("n_votus_detected")
        w.writerow(head)
        all_s = sorted(set(sample_life) | set(sample_qual) | set(abund_samples))
        for s in all_s:
            life, qual = sample_life[s], sample_qual[s]
            row = [s, sum(life.values()), life.get("lysogenic", 0),
                   life.get("lytic", 0), life.get("unknown", 0)] \
                + [qual.get(t, 0) for t in qtiers]
            if abund:
                row.append(sample_detected.get(s, 0))
            w.writerow(row)

    print(f"[build_final_table] -> {args.out_per_bin}")
    print(f"[build_final_table] -> {args.out_per_sample}")


if __name__ == "__main__":
    main()
