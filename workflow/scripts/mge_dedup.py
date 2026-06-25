#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Motor de deduplicação de elementos móveis em 2 CAMADAS + CHECKPOINT.

Problema: várias ferramentas (geNomad, IntegronFinder, MobileElementFinder,
CONJScan, ICEfinder) reportam elementos no mesmo locus. Sem dedup, o MESMO
elemento físico é contado várias vezes.

  CAMADA 1 (intra-ferramenta): junta calls redundantes DENTRO de cada tool
            (ex.: ABRicate multi-DB, CONJScan multi-modelo, hits sobrepostos).
  CAMADA 2 (entre-ferramentas): junta os calls (já limpos pela camada 1) de
            TODAS as ferramentas no mesmo locus -> 1 região não-redundante,
            anotada com quais tools a suportam (concordância = confiança).

Critério INTELIGENTE = sobreposição RECÍPROCA (cada intervalo precisa cobrir o
outro em >= overlap). Isso:
  - junta "mesmo elemento, tools diferentes" (spans parecidos);
  - NÃO junta elementos ANINHADOS (ex.: um IS dentro de um ICE), pois a cobertura
    recíproca do grande é baixa -> ficam como elementos distintos (correto).

CHECKPOINT: grava um JSON com as estatísticas das 2 camadas e uma VERIFICAÇÃO
(nenhuma sobreposição residual na saída). status=PASS só se residual==0.
"""
import argparse
import csv
import json
import os
from collections import defaultdict


# ----------------------------- geometria -----------------------------
def _ov(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def _len(iv):
    return max(1, iv[1] - iv[0] + 1)


def reciprocal_overlap(a, b):
    """min das duas frações de cobertura (0..1). Alto só se os spans são parecidos."""
    ov = _ov(a, b)
    if ov <= 0:
        return 0.0
    return min(ov / _len(a), ov / _len(b))


# ----------------------- union-find por grupo ------------------------
def _cluster(calls, overlap):
    """Agrupa calls (mesma sequência) por sobreposição recíproca. Union-find.
    n pequeno por contig -> O(n^2) é trivial e correto (independe de ordem)."""
    n = len(calls)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for i in range(n):
        ai = (calls[i]["start"], calls[i]["end"])
        for j in range(i + 1, n):
            aj = (calls[j]["start"], calls[j]["end"])
            if reciprocal_overlap(ai, aj) >= overlap:
                union(i, j)

    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(calls[i])
    return list(groups.values())


def _merge_group(members):
    start = min(m["start"] for m in members)
    end = max(m["end"] for m in members)
    tools = sorted({m["tool"] for m in members})
    types = sorted({m["element_type"] for m in members})
    score = max((float(m.get("score", 0) or 0) for m in members), default=0.0)
    return {"start": start, "end": end, "tools": tools, "types": types,
            "n_tools": len(tools), "score": score, "n_members": len(members)}


# ----------------------------- camadas -------------------------------
def layer1_intra_tool(calls, overlap):
    """Dedup DENTRO de cada (sequência, ferramenta)."""
    by_key = defaultdict(list)
    for c in calls:
        by_key[(c["seq"], c["tool"])].append(c)
    out = []
    raw_by_tool = defaultdict(int)
    l1_by_tool = defaultdict(int)
    for (seq, tool), cs in by_key.items():
        raw_by_tool[tool] += len(cs)
        for g in _cluster(cs, overlap):
            m = _merge_group(g)
            l1_by_tool[tool] += 1
            out.append({"seq": seq, "start": m["start"], "end": m["end"],
                        "tool": tool, "element_type": ";".join(m["types"]),
                        "score": m["score"]})
    return out, dict(raw_by_tool), dict(l1_by_tool)


def layer2_cross_tool(calls, overlap):
    """Junta calls (pós-camada1) de TODAS as ferramentas, por sequência."""
    by_seq = defaultdict(list)
    for c in calls:
        by_seq[c["seq"]].append(c)
    regions = []
    for seq, cs in by_seq.items():
        for g in _cluster(cs, overlap):
            m = _merge_group(g)
            m["seq"] = seq
            regions.append(m)
    return regions


def verify_no_residual(regions, overlap):
    """Auto-verificação: nenhuma região da MESMA sequência pode se sobrepor
    >= overlap (se sobrepõe, a camada 2 não deduplicou direito)."""
    by_seq = defaultdict(list)
    for r in regions:
        by_seq[r["seq"]].append(r)
    residual = 0
    for seq, rs in by_seq.items():
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                if reciprocal_overlap((rs[i]["start"], rs[i]["end"]),
                                      (rs[j]["start"], rs[j]["end"])) >= overlap:
                    residual += 1
    return residual


def run_dedup(calls, overlap):
    """Aplica as 2 camadas + verificação. Retorna (regions, checkpoint)."""
    l1, raw_by_tool, l1_by_tool = layer1_intra_tool(calls, overlap)
    regions = layer2_cross_tool(l1, overlap)
    residual = verify_no_residual(regions, overlap)
    checkpoint = {
        "overlap_threshold": overlap,
        "layer1_intra_tool": {
            "raw_calls_total": sum(raw_by_tool.values()),
            "raw_by_tool": raw_by_tool,
            "after_layer1_total": len(l1),
            "after_layer1_by_tool": l1_by_tool,
        },
        "layer2_cross_tool": {
            "input": len(l1),
            "merged_regions": len(regions),
        },
        "verification": {
            "residual_overlaps": residual,
            "status": "PASS" if residual == 0 else "FAIL",
        },
    }
    return regions, checkpoint


# --------------------------- IO / CLI --------------------------------
def read_calls(path):
    """TSV normalizado: seq, start, end, element_type, tool[, score]."""
    calls = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return calls
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                s, e = int(float(r["start"])), int(float(r["end"]))
            except (KeyError, ValueError):
                continue
            if s > e:
                s, e = e, s
            calls.append({"seq": r.get("seq", ""), "start": s, "end": e,
                          "element_type": r.get("element_type", "?"),
                          "tool": r.get("tool", "?"), "score": r.get("score", 0)})
    return calls


def write_regions(regions, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["seq", "start", "end", "n_tools", "tools",
                    "element_types", "max_score", "n_calls_merged"])
        for r in sorted(regions, key=lambda x: (x["seq"], x["start"])):
            w.writerow([r["seq"], r["start"], r["end"], r["n_tools"],
                        ";".join(r["tools"]), ";".join(r["types"]),
                        f"{r['score']:.2f}", r["n_members"]])


def _selftest():
    # 1 ICE reportado por geNomad(provirus), CONJScan, ICEfinder (~mesmo span) -> 1 regiao, 3 tools
    # 1 IS dentro do ICE (mefinder) -> NAO funde (aninhado) -> regiao separada
    # 2 hits do ABRicate no mesmo locus -> camada1 funde p/ 1
    calls = [
        {"seq": "v1", "start": 1000, "end": 18000, "element_type": "provirus", "tool": "geNomad", "score": 1},
        {"seq": "v1", "start": 1100, "end": 17800, "element_type": "conjugative_system", "tool": "CONJScan", "score": 1},
        {"seq": "v1", "start": 1050, "end": 18050, "element_type": "ICE", "tool": "ICEfinder", "score": 1},
        {"seq": "v1", "start": 8000, "end": 9200, "element_type": "IS", "tool": "mefinder", "score": 1},
        {"seq": "v1", "start": 8000, "end": 9200, "element_type": "IS", "tool": "mefinder", "score": 1},  # dup intra-tool
    ]
    regions, cp = run_dedup(calls, 0.5)
    print(json.dumps(cp, indent=2))
    print("--- regioes ---")
    for r in sorted(regions, key=lambda x: x["start"]):
        print(f"  {r['seq']} {r['start']}-{r['end']}  tools={r['tools']} types={r['types']}")
    ok = (cp["verification"]["status"] == "PASS"
          and len(regions) == 2                          # ICE(3 tools) + IS — aninhado nao fundiu
          and cp["layer1_intra_tool"]["after_layer1_by_tool"].get("mefinder") == 1)  # dup intra fundiu
    print("\nSELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", help="TSV normalizado de element calls")
    ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--out-regions")
    ap.add_argument("--out-checkpoint")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        raise SystemExit(_selftest())

    calls = read_calls(args.calls)
    regions, cp = run_dedup(calls, args.overlap)
    if args.out_regions:
        write_regions(regions, args.out_regions)
    if args.out_checkpoint:
        os.makedirs(os.path.dirname(args.out_checkpoint) or ".", exist_ok=True)
        with open(args.out_checkpoint, "w") as fh:
            json.dump(cp, fh, indent=2)
    print(f"[mge_dedup] calls={len(calls)} -> regioes={len(regions)} "
          f"| verificacao={cp['verification']['status']} "
          f"(residual={cp['verification']['residual_overlaps']})")


if __name__ == "__main__":
    main()
