#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Teste de PERMUTAÇÃO para enriquecimento de ARG/MGE na fração VIRAL vs. o
genoma completo (background).

Hipótese (Pergunta A): genes de ARG/elementos móveis estão sobre-representados
nas regiões VIRAIS (prófagos/vOTUs) em relação ao genoma inteiro?

Desenho:
  - Universo = TODOS os genes do genoma completo (background).
  - Foreground = genes nas regiões virais (is_viral=1).
  - Para cada categoria (is_arg, is_mge, ...): permutação ESTRATIFICADA por
    amostra — dentro de cada genoma, embaralha o rótulo viral/não-viral
    preservando o nº de genes virais; recomputa a contagem viral positiva;
    repete B vezes -> distribuição nula -> p empírico (enriquecimento E depleção).
  - Cruzamento paramétrico: p hipergeométrico (pooled) como sanity check.

Controla confundidor de tamanho: ao permutar RÓTULOS de genes (não posições),
o nº de genes virais é preservado, então o efeito "viral é maior" não cria
falso enriquecimento.

GUARDRAILS (anti-conclusão-silenciosa):
  - sem fração não-viral (genoma 100% viral)      -> status NO_BACKGROUND
  - poucos eventos positivos no total (< min)     -> status UNDERPOWERED
  - senão                                          -> status OK

Entrada: TSV de genes com colunas: sample, is_viral, e 1+ colunas de categoria
(is_arg, is_mge, ...). Saída: enrichment_results.tsv + JSON de checkpoint.
Stdlib apenas (roda no env pyutils).
"""
import argparse
import csv
import json
import math
import os
import random
from collections import defaultdict


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    """P(X >= k) hipergeométrico: urna N, K positivos, n sorteios."""
    lo = max(0, n - (N - K))
    hi = min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    denom = log_choose(N, n)
    total = 0.0
    for x in range(int(k), int(hi) + 1):
        total += math.exp(log_choose(K, x) + log_choose(N - K, n - x) - denom)
    return min(1.0, total)


def load_genes(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def _to01(v):
    return 1 if str(v).strip() in ("1", "True", "true", "yes", "Yes") else 0


def stratified_permutation(genes, category, B, seed=12345):
    """Permutação estratificada por sample. Retorna dict de estatísticas."""
    rng = random.Random(seed)
    by_sample = defaultdict(list)
    for g in genes:
        by_sample[g.get("sample", "ALL")].append(_to01(g.get(category, 0)))
    viral = [_to01(g.get("is_viral", 0)) for g in genes]
    by_sample_viral = defaultdict(list)
    for g in genes:
        by_sample_viral[g.get("sample", "ALL")].append(_to01(g.get("is_viral", 0)))

    # observado
    N = len(genes)
    K = sum(_to01(g.get(category, 0)) for g in genes)
    n_viral = sum(viral)
    k_obs = sum(1 for g in genes if _to01(g.get("is_viral", 0)) and _to01(g.get(category, 0)))
    n_nonviral = N - n_viral

    # nulo estratificado
    null = []
    # pré-computa, por sample, n_viral_s e o vetor de categoria
    strata = []
    for s, cats in by_sample.items():
        nv = sum(by_sample_viral[s])
        strata.append((cats, nv))
    for _ in range(B):
        k_perm = 0
        for cats, nv in strata:
            if nv == 0 or nv == len(cats):
                # nada a permutar nesse sample
                k_perm += sum(cats) if nv == len(cats) else 0
                continue
            idx = rng.sample(range(len(cats)), nv)
            k_perm += sum(cats[i] for i in idx)
        null.append(k_perm)

    ge = sum(1 for x in null if x >= k_obs)
    le = sum(1 for x in null if x <= k_obs)
    p_enrich = (1 + ge) / (B + 1)
    p_deplete = (1 + le) / (B + 1)

    exp_viral = n_viral * (K / N) if N else 0.0
    ratio = (k_obs / n_viral) / (K / N) if (n_viral and K) else float("nan")
    p_hyper = hypergeom_sf(k_obs, N, K, n_viral) if (N and K and n_viral) else float("nan")

    return {
        "category": category, "N_genes": N, "K_positive": K,
        "n_viral_genes": n_viral, "n_nonviral_genes": n_nonviral,
        "k_observed_viral_positive": k_obs,
        "expected_viral_positive": round(exp_viral, 3),
        "enrichment_ratio": (round(ratio, 3) if ratio == ratio else "NA"),
        "p_enrichment_perm": round(p_enrich, 5),
        "p_depletion_perm": round(p_deplete, 5),
        "p_enrichment_hypergeom": (round(p_hyper, 5) if p_hyper == p_hyper else "NA"),
        "permutations": B,
    }


def assess(stats, min_events):
    if stats["n_nonviral_genes"] == 0:
        return "NO_BACKGROUND"
    if stats["K_positive"] < min_events:
        return "UNDERPOWERED"
    return "OK"


def run(genes, categories, B, min_events, seed):
    results = []
    for cat in categories:
        st = stratified_permutation(genes, cat, B, seed)
        st["status"] = assess(st, min_events)
        results.append(st)
    return results


def write_results(results, out_tsv, out_json):
    cols = ["category", "status", "N_genes", "K_positive", "n_viral_genes",
            "n_nonviral_genes", "k_observed_viral_positive", "expected_viral_positive",
            "enrichment_ratio", "p_enrichment_perm", "p_depletion_perm",
            "p_enrichment_hypergeom", "permutations"]
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in results:
            w.writerow({c: r.get(c, "") for c in cols})
    if out_json:
        with open(out_json, "w") as fh:
            json.dump({"results": results,
                       "checkpoint": {"any_underpowered": any(r["status"] == "UNDERPOWERED" for r in results),
                                      "any_no_background": any(r["status"] == "NO_BACKGROUND" for r in results)}},
                      fh, indent=2)


def _selftest():
    rng = random.Random(1)
    genes = []
    # 3 amostras, 1000 genes cada, 100 virais. is_arg ENRIQUECIDO no viral
    # (10% viral vs 1% nao-viral); is_null uniforme 5% -> nao significativo.
    for s in range(3):
        for i in range(1000):
            viral = 1 if i < 100 else 0
            arg = 1 if rng.random() < (0.10 if viral else 0.01) else 0
            nul = 1 if rng.random() < 0.05 else 0
            genes.append({"sample": f"S{s}", "is_viral": viral, "is_arg": arg, "is_null": nul})
    res = run(genes, ["is_arg", "is_null"], B=2000, min_events=5, seed=7)
    for r in res:
        print(f"  {r['category']:8} status={r['status']} ratio={r['enrichment_ratio']} "
              f"p_perm={r['p_enrichment_perm']} p_hyper={r['p_enrichment_hypergeom']}")
    arg = next(r for r in res if r["category"] == "is_arg")
    nul = next(r for r in res if r["category"] == "is_null")
    ok = (arg["status"] == "OK" and arg["p_enrichment_perm"] < 0.01
          and float(arg["enrichment_ratio"]) > 3
          and nul["p_enrichment_perm"] > 0.05)
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", help="TSV: sample, is_viral, is_arg, is_mge, ...")
    ap.add_argument("--categories", help="virgula-separado; default = colunas is_* (exceto is_viral)")
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--min-events", type=int, default=5)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out-tsv")
    ap.add_argument("--out-json")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        raise SystemExit(_selftest())

    genes = load_genes(args.genes)
    if args.categories:
        cats = args.categories.split(",")
    else:
        sample_keys = genes[0].keys() if genes else []
        cats = [c for c in sample_keys if c.startswith("is_") and c != "is_viral"]
    results = run(genes, cats, args.permutations, args.min_events, args.seed)
    write_results(results, args.out_tsv or "enrichment_results.tsv", args.out_json)
    for r in results:
        print(f"[enrichment] {r['category']}: status={r['status']} "
              f"ratio={r['enrichment_ratio']} p_perm={r['p_enrichment_perm']}")


if __name__ == "__main__":
    main()
