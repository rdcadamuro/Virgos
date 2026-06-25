#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Consolida TODOS os sinais de mobiloma + DEDUP de 2 CAMADAS por coordenada.

Duas saidas complementares:

(1) TABELA POR vOTU  -> mobilome_table.tsv (1 linha por vOTU; contadores +
    'mobilome_flag' + 'n_distinct_mge_regions'). Visao "quem tem o que".

(2) REGIOES NAO-REDUNDANTES -> mobile_element_regions.tsv + checkpoint JSON.
    Cada ferramenta emite ELEMENT CALLS com coordenada no FRAME do vOTU:
       seq, start, end, element_type, tool, score
    Fontes de calls:
       - geNomad   : provirus (vOTU inteiro) + integrase/conjugacao/transposase
                     com coordenada (build_mobilome.py --out-calls)
       - IntegronFinder : integrons (.integrons: pos_beg/pos_end -> span por integron)
       - MobileElementFinder : IS/transposons (CSV: start/end)
       - CONJScan  : sistemas conjugativos (run_conjscan.py, span por sistema)
       - ICEfinder (opcional): ICEs, se o gancho gravar calls
    Todas as calls passam pelo motor mge_dedup (2 camadas + verificacao). O
    checkpoint so e PASS se nao sobra sobreposicao residual.

tRNAscan entra como MARCADOR (n_tRNAs na tabela), NAO como elemento.
Parsers tolerantes: arquivo ausente/vazio = 0 (degradacao graciosa).
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mge_dedup import run_dedup, write_regions  # motor validado


def read_fasta_names(path):
    out = []
    if os.path.exists(path):
        with open(path, errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(">"):
                    out.append(ln[1:].split()[0].strip())
    return out


def load_genomad(path):
    rows = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path, errors="ignore", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for d in r:
            rows[d.get("votu_rep", "")] = d
    return rows


def read_calls_tsv(path):
    """Le um TSV de element calls (seq,start,end,element_type,tool,score)."""
    out = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path, errors="ignore", newline="") as fh:
        for d in csv.DictReader(fh, delimiter="\t"):
            try:
                s, e = int(float(d["start"])), int(float(d["end"]))
            except (KeyError, ValueError):
                continue
            out.append({"seq": d.get("seq", ""), "start": s, "end": e,
                        "element_type": d.get("element_type", "?"),
                        "tool": d.get("tool", "?"), "score": d.get("score", 0)})
    return out


# ---------- parsers que tambem geram contagem por-vOTU + calls --------------
def parse_integron(path):
    """IntegronFinder .summary -> {contig: (n_total, n_complete)} (contagem)."""
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path, errors="ignore") as fh:
        header, idx = None, {}
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"):
                continue
            f = ln.split("\t")
            if header is None and "ID_replicon" in f:
                header = f
                idx = {c: i for i, c in enumerate(f)}
                continue
            if header is None:
                continue
            contig = f[idx.get("ID_replicon", 0)]

            def gi(col):
                try:
                    return int(float(f[idx[col]])) if col in idx and f[idx[col]] not in ("", "NA") else 0
                except (ValueError, IndexError):
                    return 0
            comp, in0, calin = gi("complete"), gi("In0"), gi("CALIN")
            out[contig] = (comp + in0 + calin, comp)
    return out


def integron_calls(integron_dir):
    """Varre os .integrons do IntegronFinder -> 1 call por integron (span dos
    elementos), no frame do vOTU (IntegronFinder roda no votus.fna)."""
    calls = []
    if not integron_dir or not os.path.isdir(integron_dir):
        return calls
    for root, _dirs, files in os.walk(integron_dir):
        for fn in files:
            if not fn.endswith(".integrons"):
                continue
            fp = os.path.join(root, fn)
            try:
                with open(fp, errors="ignore") as fh:
                    header, idx = None, {}
                    spans = {}  # (replicon, ID_integron) -> [pos...]
                    for ln in fh:
                        ln = ln.rstrip("\n")
                        if not ln or ln.startswith("#"):
                            continue
                        f = ln.split("\t")
                        if header is None and ("ID_replicon" in f or "pos_beg" in f):
                            header = f
                            idx = {c: i for i, c in enumerate(f)}
                            continue
                        if header is None:
                            continue

                        def gv(col):
                            i = idx.get(col)
                            return f[i] if i is not None and i < len(f) else ""
                        rep = gv("ID_replicon")
                        iid = gv("ID_integron") or "1"
                        try:
                            pb = int(float(gv("pos_beg"))); pe = int(float(gv("pos_end")))
                        except ValueError:
                            continue
                        spans.setdefault((rep, iid), []).extend([pb, pe])
                    for (rep, _iid), pts in spans.items():
                        if rep and pts:
                            calls.append({"seq": rep, "start": min(pts), "end": max(pts),
                                          "element_type": "integron", "tool": "IntegronFinder",
                                          "score": 1})
            except OSError:
                continue
    return calls


def parse_mefinder(path):
    """MobileElementFinder .csv -> ({contig:(count,types)} , [calls]).
    Calls usam as colunas de coordenada do MGEdb (start/end no contig)."""
    counts, calls = {}, []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return counts, calls
    with open(path, errors="ignore", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return counts, calls
    rows = list(csv.reader(lines))
    if not rows:
        return counts, calls
    header = rows[0]
    low = [h.strip().lower() for h in header]

    def col(*names):
        for nm in names:
            if nm in low:
                return low.index(nm)
        return None
    ci = col("contig")
    if ci is None:
        ci = next((i for i, h in enumerate(low) if "contig" in h or "sequence" in h), None)
    ni = col("name", "type", "mge")
    si = col("start", "contig_start", "begin")
    ei = col("end", "contig_end", "stop")
    if ci is None:
        return counts, calls
    for r in rows[1:]:
        if len(r) <= ci or not r[ci].strip():
            continue
        contig = r[ci].strip()
        cnt, types = counts.get(contig, (0, set()))
        cnt += 1
        if ni is not None and len(r) > ni and r[ni].strip():
            types.add(r[ni].strip())
        counts[contig] = (cnt, types)
        if si is not None and ei is not None and len(r) > max(si, ei):
            try:
                s, e = int(float(r[si])), int(float(r[ei]))
                calls.append({"seq": contig, "start": s, "end": e,
                              "element_type": "IS_transposon", "tool": "mefinder",
                              "score": 1})
            except ValueError:
                pass
    return counts, calls


def parse_trnascan(path):
    """tRNAscan-SE -o -> {contig: n_tRNAs} (so MARCADOR; nao vira elemento)."""
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    started = False
    with open(path, errors="ignore") as fh:
        for ln in fh:
            if ln.startswith("----") or ln.startswith("--------"):
                started = True
                continue
            if not started:
                if ln.startswith(("Sequence", "Name")) or not ln.strip():
                    continue
            f = ln.rstrip("\n").split("\t")
            if not f or not f[0].strip():
                continue
            contig = f[0].strip()
            out[contig] = out.get(contig, 0) + 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--votus", required=True)
    ap.add_argument("--genomad", required=True, help="mobilome_genomad.tsv (tabela por vOTU)")
    ap.add_argument("--integron", required=True, help=".summary do IntegronFinder")
    ap.add_argument("--mefinder", required=True)
    ap.add_argument("--trnascan", required=True)
    ap.add_argument("--icefinder", required=True)
    ap.add_argument("--out", required=True, help="mobilome_table.tsv (por vOTU)")
    # element calls / dedup
    ap.add_argument("--genomad-calls", default="", help="calls do build_mobilome --out-calls")
    ap.add_argument("--conjscan-calls", default="", help="calls do run_conjscan.py")
    ap.add_argument("--icefinder-calls", default="", help="(opcional) calls do ICEfinder")
    ap.add_argument("--integron-dir", default="", help="dir com os .integrons (coords)")
    ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--out-regions", default="", help="mobile_element_regions.tsv")
    ap.add_argument("--out-checkpoint", default="", help="mobilome_dedup_checkpoint.json")
    args = ap.parse_args()

    reps = read_fasta_names(args.votus)
    gen = load_genomad(args.genomad)
    integ = parse_integron(args.integron)
    mef_counts, mef_calls = parse_mefinder(args.mefinder)
    trna = parse_trnascan(args.trnascan)

    # ----- ELEMENT CALLS de todas as ferramentas (frame do vOTU) -----
    all_calls = []
    all_calls += read_calls_tsv(args.genomad_calls)
    all_calls += read_calls_tsv(args.conjscan_calls)
    all_calls += read_calls_tsv(args.icefinder_calls)
    all_calls += integron_calls(args.integron_dir)
    all_calls += mef_calls

    # ----- DEDUP 2 camadas + checkpoint -----
    regions, checkpoint = run_dedup(all_calls, args.overlap)
    n_reg_by_seq = {}
    for r in regions:
        n_reg_by_seq[r["seq"]] = n_reg_by_seq.get(r["seq"], 0) + 1
    if args.out_regions:
        write_regions(regions, args.out_regions)
    if args.out_checkpoint:
        import json
        os.makedirs(os.path.dirname(args.out_checkpoint) or ".", exist_ok=True)
        with open(args.out_checkpoint, "w") as fh:
            json.dump(checkpoint, fh, indent=2)

    # ----- TABELA por vOTU -----
    cols = ["votu_rep", "sample", "topology", "is_provirus",
            "n_amr_genes", "amr_genes", "n_conjugation_genes",
            "n_integrase", "n_transposase",
            "n_integrons", "n_complete_integrons",
            "n_mobile_elements", "mobile_element_types",
            "n_tRNAs", "n_distinct_mge_regions", "mobilome_flag"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for rep in reps:
            g = gen.get(rep, {})
            n_int, n_intc = integ.get(rep, (0, 0))
            n_mef, mef_types = mef_counts.get(rep, (0, set()))
            n_trna = trna.get(rep, 0)
            n_reg = n_reg_by_seq.get(rep, 0)
            is_prov = (g.get("is_provirus", "No") == "Yes")
            n_conj = int(g.get("n_conjugation_genes", 0) or 0)
            n_integrase = int(g.get("n_integrase", 0) or 0)
            flag = any([is_prov, n_int > 0, n_mef > 0, n_conj > 0, n_integrase > 0, n_reg > 0])
            w.writerow([
                rep, g.get("sample", ""), g.get("topology", ""),
                "Yes" if is_prov else "No",
                g.get("n_amr_genes", 0), g.get("amr_genes", "-"), n_conj,
                n_integrase, g.get("n_transposase", 0),
                n_int, n_intc,
                n_mef, ";".join(sorted(mef_types)) if mef_types else "-",
                n_trna, n_reg, "Yes" if flag else "No",
            ])

    ver = checkpoint.get("verification", {})
    print(f"[consolidate_mobilome] {len(reps)} vOTUs -> {args.out} | "
          f"calls={len(all_calls)} -> regioes={len(regions)} "
          f"(dedup {ver.get('status')}, residual={ver.get('residual_overlaps')})")


if __name__ == "__main__":
    main()
