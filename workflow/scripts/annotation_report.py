#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relatorio consolidado da anotacao funcional:

  (1) TSV GERAL detalhado por vOTU  -> annotation_master.tsv
      junta mobiloma (por vOTU) + nº de loci de resistoma (ABRicate) + nº de
      regioes MGE nao-redundantes (apos dedup 2 camadas).

  (2) LOG ESTATISTICO legivel        -> annotation_stats.log
      - enriquecimento POR CATEGORIA (ARG, conjugacao, integrase, transposase,
        integron, ICE/conjugative_system, IS, MGE) com ratio, p (permutacao +
        hipergeometrica), status e INTERPRETACAO automatica;
      - checkpoint da dedup (2 camadas) — PASS/FAIL + residual;
      - totais de resistoma e de mobiloma.

Defensivo: qualquer entrada ausente/vazia -> secao marcada como indisponivel.
"""
import argparse
import csv
import json
import os


def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, errors="ignore", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return {}


def interpret(r):
    """Interpretacao automatica de uma linha de enriquecimento."""
    st = r.get("status", "?")
    if st == "NO_BACKGROUND":
        return "sem background nao-viral (genoma 100% viral) — N/A"
    if st == "UNDERPOWERED":
        return f"poucos eventos (K={r.get('K_positive')}) — sem poder estatistico"
    try:
        pe = float(r.get("p_enrichment_perm", 1))
        pd = float(r.get("p_depletion_perm", 1))
        ratio = float(r.get("enrichment_ratio")) if r.get("enrichment_ratio") not in ("", "NA") else float("nan")
    except (TypeError, ValueError):
        return "indeterminado"
    if pe <= 0.05:
        return f"ENRIQUECIDO no viral (ratio={ratio:.2f}, p={pe:.4g})"
    if pd <= 0.05:
        return f"DEPLETADO no viral (ratio={ratio:.2f}, p={pd:.4g})"
    return f"sem enriquecimento (ratio={ratio:.2f}, p={pe:.3g})"


def build_master(annot_dir, out_tsv):
    mob = read_tsv(os.path.join(annot_dir, "mobilome_table.tsv"))
    res = read_tsv(os.path.join(annot_dir, "resistome_abricate.tsv"))
    regions = read_tsv(os.path.join(annot_dir, "mobile_element_regions.tsv"))

    # loci de resistoma por vOTU (coluna 'sequence')
    res_by_votu = {}
    res_genes = {}
    for r in res:
        sq = r.get("sequence", "")
        res_by_votu[sq] = res_by_votu.get(sq, 0) + 1
        g = r.get("best_gene", "") or r.get("all_genes", "")
        if g:
            res_genes.setdefault(sq, []).append(g)
    # tools concordantes por vOTU (das regioes)
    reg_tools = {}
    for r in regions:
        sq = r.get("seq", "")
        reg_tools.setdefault(sq, set()).update(
            (r.get("tools", "") or "").split(";") if r.get("tools") else [])

    base_cols = (list(mob[0].keys()) if mob else
                 ["votu_rep", "sample", "mobilome_flag"])
    extra = ["n_arg_loci_abricate", "arg_genes_abricate", "mge_region_tools"]
    cols = base_cols + extra
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for m in mob:
            rep = m.get("votu_rep", "")
            row = [m.get(c, "") for c in base_cols]
            row += [res_by_votu.get(rep, 0),
                    ";".join(sorted(set(res_genes.get(rep, [])))) or "-",
                    ";".join(sorted(t for t in reg_tools.get(rep, set()) if t)) or "-"]
            w.writerow(row)
    return len(mob), sum(res_by_votu.values()), len(regions)


def build_stats_log(annot_dir, out_log, n_votus, n_res, n_regions):
    enr = read_tsv(os.path.join(annot_dir, "enrichment_results.tsv"))
    dedup_cp = load_json(os.path.join(annot_dir, "mobilome_dedup_checkpoint.json"))
    mob = read_tsv(os.path.join(annot_dir, "mobilome_table.tsv"))

    L = []
    L.append("=" * 70)
    L.append("  RELATORIO ESTATISTICO — ANOTACAO FUNCIONAL (resistoma + mobiloma)")
    L.append(f"  diretorio: {annot_dir}")
    L.append("=" * 70)

    # --- enriquecimento por categoria ---
    L.append("")
    L.append("[ENRIQUECIMENTO — ARG/MGE na fracao VIRAL vs GENOMA COMPLETO]")
    if not enr:
        L.append("  (sem resultados — genoma sem genes anotados ou teste nao rodado)")
    else:
        e0 = enr[0]
        L.append(f"  universo: {e0.get('N_genes')} genes | viral: {e0.get('n_viral_genes')} | "
                 f"nao-viral: {e0.get('n_nonviral_genes')} | permutacoes: {e0.get('permutations')}")
        L.append("  " + "-" * 66)
        L.append(f"  {'categoria':<22}{'obs/esp':<12}{'ratio':<8}{'p_perm':<9}{'status':<13}interpretacao")
        for r in enr:
            cat = r.get("category", "?")
            obs = f"{r.get('k_observed_viral_positive')}/{r.get('expected_viral_positive')}"
            ratio = r.get("enrichment_ratio", "NA")
            pe = r.get("p_enrichment_perm", "NA")
            st = r.get("status", "?")
            L.append(f"  {cat:<22}{obs:<12}{str(ratio):<8}{str(pe):<9}{st:<13}{interpret(r)}")

    # --- dedup ---
    L.append("")
    L.append("[DEDUP DE MOBILOMA — 2 camadas (intra-tool + entre-tools)]")
    if dedup_cp:
        l1 = dedup_cp.get("layer1_intra_tool", {})
        l2 = dedup_cp.get("layer2_cross_tool", {})
        ver = dedup_cp.get("verification", {})
        L.append(f"  camada1: {l1.get('raw_calls_total', '?')} calls -> {l1.get('after_layer1_total', '?')} "
                 f"(por tool: {l1.get('raw_by_tool', {})})")
        L.append(f"  camada2: {l2.get('input', '?')} -> {l2.get('merged_regions', '?')} regioes nao-redundantes")
        L.append(f"  VERIFICACAO: status={ver.get('status')} residual_overlaps={ver.get('residual_overlaps')}")
        if ver.get("status") != "PASS":
            L.append("  *** ATENCAO: dedup NAO verificada (sobreposicao residual) — investigar ***")
    else:
        L.append("  (sem checkpoint de dedup)")

    # --- totais ---
    L.append("")
    L.append("[TOTAIS]")
    n_flag = sum(1 for m in mob if m.get("mobilome_flag") == "Yes")
    L.append(f"  vOTUs: {n_votus} | com sinal de elemento movel (mobilome_flag=Yes): {n_flag}")
    L.append(f"  loci de ARG/virulencia (ABRicate, dedup): {n_res}")
    L.append(f"  regioes MGE nao-redundantes (apos dedup): {n_regions}")
    n_prov = sum(1 for m in mob if m.get("is_provirus") == "Yes")
    L.append(f"  vOTUs provirus: {n_prov}")
    L.append("=" * 70)

    os.makedirs(os.path.dirname(out_log) or ".", exist_ok=True)
    text = "\n".join(L) + "\n"
    with open(out_log, "w") as fh:
        fh.write(text)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot-dir", required=True, help="diretorio 0X_annotation")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-log", required=True)
    args = ap.parse_args()

    n_votus, n_res, n_regions = build_master(args.annot_dir, args.out_tsv)
    text = build_stats_log(args.annot_dir, args.out_log, n_votus, n_res, n_regions)
    print(text)
    print(f"[annotation_report] -> {args.out_tsv} | {args.out_log}")


if __name__ == "__main__":
    main()
