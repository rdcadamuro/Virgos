#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verificacao INTELIGENTE de execucao do pipeline (anti-erro-silencioso).

Roda no fim do pipeline (ou a qualquer momento) e produz um relatorio claro:
  - confirma que CADA etapa COMPLETOU (le os marcadores "OK -- rodou ate o fim");
  - conta os resultados de cada etapa (virus, contigs, bins, vOTUs, taxonomia...);
  - CRUZA a consistencia entre etapas para achar falhas SILENCIOSAS
    (ex.: geNomad achou 50 virus mas o CheckV avaliou 0 -> algo quebrou no meio);
  - VARRE os logs por assinaturas de erro que possam ter passado batido
    (Traceback, "No module named", "command not found", core dumped, OOM...).

Classifica cada item como:  OK | VAZIO-OK (vazio legitimo) | AVISO | FALHA.
Veredito final: INTEGRO  /  INTEGRO-COM-AVISOS  /  FALHA(S) DETECTADA(S).

Uso:
  python verify_run.py --outdir <saida> --mode genome|metagenome \
      --samples <samples.tsv> --out <saida>/.../RUN_HEALTHCHECK.txt [--strict]
"""
import argparse
import json
import os
import sys

# Assinaturas de erro de ALTA PRECISAO (raramente aparecem fora de falha real)
FATAL_SIGNATURES = [
    "Traceback (most recent call last)",
    "*** ERRO REAL",
    "No module named",
    "command not found",
    "CalledProcessError",
    "Segmentation fault",
    "core dumped",
    "MemoryError",
    "Cannot allocate",
    "cannot allocate",
    "Killed",
]

OK, EMPTY_OK, WARN, FAIL = "OK", "VAZIO-OK", "AVISO", "FALHA"
SYM = {OK: "[ OK ]", EMPTY_OK: "[VAZIO]", WARN: "[AVISO]", FAIL: "[FALHA]"}


def count_fasta(path):
    if not path or not os.path.exists(path):
        return None
    n = 0
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def count_rows(path):
    """Linhas de dados (sem header, sem linhas em branco)."""
    if not path or not os.path.exists(path):
        return None
    n = 0
    with open(path, errors="ignore") as fh:
        next(fh, None)  # header
        for line in fh:
            if line.strip():
                n += 1
    return n


def scan_log(path):
    """Retorna lista de (linha, texto) com assinaturas fatais no log."""
    hits = []
    if not path or not os.path.exists(path):
        return hits
    with open(path, errors="ignore") as fh:
        for i, line in enumerate(fh, 1):
            for sig in FATAL_SIGNATURES:
                if sig in line:
                    hits.append((i, line.rstrip()[:200]))
                    break
    return hits


def read_samples(path):
    """Retorna lista de dicts: sample, read_type, has_reads."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0]:
                continue
            rt = p[idx["read_type"]] if "read_type" in idx and len(p) > idx["read_type"] else ""
            r1 = p[idx["reads1"]] if "reads1" in idx and len(p) > idx["reads1"] else ""
            rows.append({
                "sample": p[0],
                "read_type": rt,
                "has_reads": rt in {"ilu-paired", "ilu-single", "nano"} and bool(r1),
            })
    return rows


class Report:
    def __init__(self):
        self.lines = []
        self.counts = {OK: 0, EMPTY_OK: 0, WARN: 0, FAIL: 0}

    def item(self, status, stage, msg):
        self.counts[status] += 1
        self.lines.append(f"  {SYM[status]} {stage:<16} {msg}")

    def section(self, title):
        self.lines.append("")
        self.lines.append(title)

    def raw(self, text):
        self.lines.append(text)


def check_log(rep, stage, log_path, completed_marker=None):
    """Varre o log por erros fatais e (opcional) confirma marcador de termino.
    Retorna True se achou erro fatal."""
    hits = scan_log(log_path)
    if hits:
        ln, txt = hits[0]
        rep.item(FAIL, stage, f"log com sinal de ERRO (linha {ln}): {txt}")
        return True
    if completed_marker and os.path.exists(log_path):
        with open(log_path, errors="ignore") as fh:
            if completed_marker not in fh.read():
                rep.item(WARN, stage,
                         f"log sem o marcador de termino ('{completed_marker}') — verifique")
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mode", choices=["genome", "metagenome"], required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true",
                    help="sai com codigo !=0 se houver FALHA (aborta o pipeline)")
    args = ap.parse_args()

    o = args.outdir
    D_GENOMAD = f"{o}/01_genomad"
    D_CHECKV = f"{o}/02_checkv"
    D_BINS = f"{o}/03_bins_vrhyme_vamb"
    D_PHABOX = f"{o}/04_phabox2"
    D_VITAP = f"{o}/05_vitap"
    D_ABUND = f"{o}/06_votu_abundance"
    if args.mode == "metagenome":
        D_ANNOT = f"{o}/07_annotation"
        D_SUM = f"{o}/08_summary_final"
    else:
        D_ANNOT = f"{o}/06_annotation"
        D_SUM = f"{o}/07_summary_final"
    LOGS = f"{o}/logs"

    samples = read_samples(args.samples)
    rep = Report()
    rep.raw("=" * 64)
    rep.raw("  VERIFICACAO INTELIGENTE DE EXECUCAO")
    rep.raw(f"  saida: {o}")
    rep.raw(f"  modo : {args.mode}   |   amostras: {len(samples)}")
    rep.raw("=" * 64)

    # ---------------- por amostra: geNomad + CheckV ----------------
    for s in samples:
        sn = s["sample"]
        rep.section(f"[amostra: {sn}]  (read_type={s['read_type'] or 'sem reads'})")

        # 01 geNomad
        viral = f"{D_GENOMAD}/viral_contigs/{sn}.viral.fna"
        n_viral = count_fasta(viral)
        glog = f"{LOGS}/01_genomad/{sn}.log"
        gfail = check_log(rep, "01 geNomad", glog, "rodou ate o fim")
        if n_viral is None:
            rep.item(FAIL, "01 geNomad", f"saida AUSENTE: {viral}")
            n_viral = 0
        elif n_viral == 0:
            rep.item(EMPTY_OK, "01 geNomad", "0 contigs virais (geNomad completou; amostra sem virus)")
        else:
            rep.item(OK, "01 geNomad", f"{n_viral} contigs virais identificados")

        # 02 CheckV — caminho depende de ter reads (bins) ou nao (contigs)
        if s["has_reads"]:
            csum = f"{D_CHECKV}/{sn}_bins/quality_summary.tsv"
            clog = f"{LOGS}/02_checkv_bins/{sn}.log"
            label = "02 CheckV-bins"
            linked = f"{D_BINS}/{sn}/linked_bins.fna"
            n_units = count_fasta(linked)  # nº de bins linkados esperados
        else:
            csum = f"{D_CHECKV}/{sn}/quality_summary.tsv"
            clog = f"{LOGS}/02_checkv/{sn}.log"
            label = "02 CheckV"
            n_units = n_viral  # CheckV avalia os contigs virais

        cfail = check_log(rep, label, clog, "rodou ate o fim")
        n_check = count_rows(csum)
        if n_check is None:
            rep.item(FAIL, label, f"saida AUSENTE: {csum}")
        elif n_check == 0:
            if n_units and not cfail:
                # *** padrao classico de FALHA SILENCIOSA ***
                rep.item(FAIL, label,
                         f"INCONSISTENCIA: havia {n_units} unidades para avaliar, "
                         f"mas o CheckV produziu 0 linhas. Provavel falha silenciosa.")
            elif n_units and cfail:
                pass  # erro ja reportado pela varredura do log; nao duplicar
            else:
                rep.item(EMPTY_OK, label, "0 avaliacoes (nada a avaliar a montante)")
        else:
            rep.item(OK, label, f"{n_check} unidades avaliadas (completude/qualidade)")

        # 03 binning (so com reads): confere vRhyme e VAMB
        if s["has_reads"]:
            for binner, blog in (("vRhyme", f"{LOGS}/03_vrhyme/{sn}.log"),
                                  ("VAMB", f"{LOGS}/03_vamb/{sn}.log")):
                bf = check_log(rep, f"03 {binner}", blog, "rodou ate o fim")
                if not bf and os.path.exists(blog):
                    with open(blog, errors="ignore") as fh:
                        txt = fh.read()
                    if "rodou ate o fim" in txt or "0 bins" in txt or "NENHUM bin" in txt:
                        # extrai a linha de status do nosso echo
                        status = next((l for l in txt.splitlines()
                                       if l.startswith(f"[{binner}]")), "")
                        rep.item(OK, f"03 {binner}", status.strip() or "rodou ate o fim")

    # ---------------- globais: derep / phabox / vitap / abund / summary --------
    rep.section("[etapas globais]")

    votus = f"{D_BINS}/derep/votus.fna"
    n_votus = count_fasta(votus)
    dlog = f"{LOGS}/03_derep.log"
    check_log(rep, "03 derep", dlog)
    if n_votus is None:
        rep.item(FAIL, "03 derep", f"saida AUSENTE: {votus}")
        n_votus = 0
    elif n_votus == 0:
        rep.item(EMPTY_OK, "03 derep", "0 vOTUs (nenhuma unidade viral aprovada a montante)")
    else:
        rep.item(OK, "03 derep", f"{n_votus} vOTUs (representantes dereplicados)")

    # 04 PhaBox2
    phatyp = f"{D_PHABOX}/final_prediction/phatyp_prediction.tsv"
    plog = f"{LOGS}/04_phabox.log"
    pfail = check_log(rep, "04 PhaBox2", plog, "rodou ate o fim")
    n_pht = count_rows(phatyp)
    if n_pht is None:
        rep.item(FAIL, "04 PhaBox2", f"saida AUSENTE: {phatyp}")
    elif n_pht == 0:
        if n_votus > 0 and not pfail:
            rep.item(FAIL, "04 PhaBox2",
                     f"INCONSISTENCIA: {n_votus} vOTUs mas 0 predicoes de lifestyle. Falha silenciosa?")
        elif n_votus > 0 and pfail:
            pass  # erro ja reportado pela varredura do log
        else:
            rep.item(EMPTY_OK, "04 PhaBox2", "0 predicoes (sem vOTUs a montante)")
    else:
        rep.item(OK, "04 PhaBox2", f"{n_pht} predicoes de lifestyle")

    # 05 VITAP
    vtax = f"{D_VITAP}/vitap_taxonomy.tsv"
    vlog = f"{LOGS}/05_vitap.log"
    vfail = check_log(rep, "05 VITAP", vlog, "rodou ate o fim")
    n_vt = count_rows(vtax)
    if n_vt is None:
        rep.item(FAIL, "05 VITAP", f"saida AUSENTE: {vtax}")
    elif n_votus > 0 and n_vt == 0:
        # VITAP pode legitimamente dar 0 p/ fragmentos divergentes -> AVISO, nao FALHA
        rep.item(WARN, "05 VITAP",
                 f"{n_votus} vOTUs mas 0 linhagens. Pode ser legitimo (fragmentos "
                 f"curtos/divergentes do ICTV); confirme no log se VITAP completou.")
    elif n_vt == 0:
        rep.item(EMPTY_OK, "05 VITAP", "0 linhagens (sem vOTUs a montante)")
    else:
        rep.item(OK, "05 VITAP", f"{n_vt} linhagens atribuidas")

    # 06 abundancia (so metagenoma)
    if args.mode == "metagenome":
        tpm = f"{D_ABUND}/votu_abundance_tpm.tsv"
        n_ab = count_rows(tpm)
        if n_ab is None:
            rep.item(FAIL, "06 abundancia", f"saida AUSENTE: {tpm}")
        elif n_votus > 0 and n_ab == 0:
            rep.item(FAIL, "06 abundancia",
                     f"INCONSISTENCIA: {n_votus} vOTUs mas matriz de abundancia vazia.")
        elif n_ab == 0:
            rep.item(EMPTY_OK, "06 abundancia", "matriz vazia (sem vOTUs)")
        else:
            rep.item(OK, "06 abundancia", f"{n_ab} vOTUs na matriz de abundancia")
        if any(s["read_type"] == "nano" for s in samples):
            rep.raw("  [nota]  Nanopore detectado -> limiar de identidade da abundancia "
                    "relaxado p/ nano (abundance.min_read_pct_id_nano). Reads long-read "
                    "tem mais erro; 95% de Illumina subestimaria a abundancia.")

    # 06/07 anotacao: resistoma (ABRicate) + mobiloma (geNomad)
    res = f"{D_ANNOT}/resistome_abricate.tsv"
    mob = f"{D_ANNOT}/mobilome_table.tsv"
    check_log(rep, "06 resistoma", f"{LOGS}/06_abricate.log", "rodou ate o fim")
    n_res = count_rows(res)
    by_db = f"{D_ANNOT}/resistome_by_db.tsv"
    n_db = count_rows(by_db)  # nº de bancos comparados (do resumo por banco)
    if n_res is None:
        rep.item(FAIL, "06 resistoma", f"saida AUSENTE: {res}")
    else:
        extra = f" | {n_db} bancos comparados" if n_db else ""
        rep.item(OK, "06 resistoma",
                 f"{n_res} loci de ARG/virulencia/plasmideo (ABRicate dedup){extra}. 0 e comum em fagos.")
    # 06 resistoma: checkpoint de DEDUP por locus — verifica residual==0
    res_cp = f"{D_ANNOT}/resistome_dedup_checkpoint.json"
    if os.path.exists(res_cp):
        try:
            with open(res_cp) as fh:
                rj = json.load(fh)
            ver = rj.get("verification", {})
            if ver.get("status") == "PASS" and ver.get("residual_overlaps", 1) == 0:
                rep.item(OK, "06 dedup-ARG",
                         f"loci nao-redundantes; {rj.get('n_raw_hits', '?')} hits -> "
                         f"{rj.get('n_loci', '?')} loci, residual=0")
            else:
                rep.item(FAIL, "06 dedup-ARG",
                         f"DEDUP NAO VERIFICADA: residual={ver.get('residual_overlaps')} status={ver.get('status')}")
        except (ValueError, OSError) as e:
            rep.item(WARN, "06 dedup-ARG", f"checkpoint ilegivel: {e}")
    cmob = check_log(rep, "06 mobiloma", f"{LOGS}/06_mobilome.log", "rodou ate o fim")
    n_mob = count_rows(mob)
    if n_mob is None:
        rep.item(FAIL, "06 mobiloma", f"saida AUSENTE: {mob}")
    elif n_votus > 0 and n_mob == 0 and not cmob:
        rep.item(FAIL, "06 mobiloma", f"INCONSISTENCIA: {n_votus} vOTUs mas mobiloma vazio.")
    elif n_mob > 0:
        rep.item(OK, "06 mobiloma", f"{n_mob} vOTUs anotados (provirus/AMR/MGE)")

    # 06 enriquecimento (permutacao ARG/MGE viral vs genoma) — le o checkpoint
    enr_cp = f"{D_ANNOT}/enrichment_checkpoint.json"
    if os.path.exists(enr_cp):
        try:
            with open(enr_cp) as fh:
                ej = json.load(fh)
            results = ej.get("results", [])
            if not results:
                rep.item(EMPTY_OK, "06 enriquecimento",
                         "sem teste (genoma sem genes anotados ou sem background nao-viral)")
            for r in results:
                st = r.get("status", "?")
                msg = (f"{r.get('category')}: ratio={r.get('enrichment_ratio')} "
                       f"p_perm={r.get('p_enrichment_perm')} [{st}]")
                if st == "OK":
                    rep.item(OK, "06 enriquecimento", msg)
                else:
                    rep.item(WARN, "06 enriquecimento",
                             msg + " — NO_BACKGROUND/UNDERPOWERED: interpretar com cautela")
        except (ValueError, OSError) as e:
            rep.item(WARN, "06 enriquecimento", f"checkpoint ilegivel: {e}")

    # 06 mobiloma: checkpoint de DEDUP (2 camadas) — verifica residual==0
    dedup_cp = f"{D_ANNOT}/mobilome_dedup_checkpoint.json"
    if os.path.exists(dedup_cp):
        try:
            with open(dedup_cp) as fh:
                dj = json.load(fh)
            ver = dj.get("verification", {})
            if ver.get("status") == "PASS" and ver.get("residual_overlaps", 1) == 0:
                rep.item(OK, "06 dedup-MGE",
                         f"2 camadas OK; regioes={dj.get('layer2_cross_tool', {}).get('merged_regions', '?')}, residual=0")
            else:
                rep.item(FAIL, "06 dedup-MGE",
                         f"DEDUP NAO VERIFICADA: residual={ver.get('residual_overlaps')} status={ver.get('status')}")
        except (ValueError, OSError) as e:
            rep.item(WARN, "06 dedup-MGE", f"checkpoint ilegivel: {e}")

    # summary final
    per_bin = f"{D_SUM}/viral_bins_table.tsv"
    per_sample = f"{D_SUM}/per_sample_summary.tsv"
    for name, path in (("tabela vOTUs", per_bin), ("resumo amostras", per_sample)):
        n = count_rows(path)
        if n is None:
            rep.item(FAIL, "07 summary", f"saida AUSENTE: {path}")
        else:
            rep.item(OK, "07 summary", f"{name}: {n} linhas")

    # 09 ENA: pacote de submissao MIUViG (se gerado). Confere consistencia
    # uViGs<->manifestos e conta os campos TODO que o submissor precisa preencher.
    ena_meta = f"{D_SUM}/ena_submission/ena_uvig_metadata.tsv"
    if os.path.exists(ena_meta):
        n_uvig = count_rows(ena_meta)
        manif_dir = f"{D_SUM}/ena_submission/manifests"
        n_manif = len([f for f in os.listdir(manif_dir) if f.endswith(".manifest")]) \
            if os.path.isdir(manif_dir) else 0
        n_todo = 0
        try:
            with open(ena_meta, errors="ignore") as fh:
                n_todo = sum(ln.count("TODO_PREENCHER") for ln in fh)
        except OSError:
            pass
        if n_uvig and n_manif == n_uvig:
            rep.item(OK, "09 ENA",
                     f"{n_uvig} uViGs = {n_manif} manifestos MIUViG; "
                     f"{n_todo} campos TODO a preencher (data/local/extracao/assembler)")
        elif n_uvig:
            rep.item(WARN, "09 ENA",
                     f"INCONSISTENCIA: {n_uvig} uViGs mas {n_manif} manifestos")
        else:
            rep.item(EMPTY_OK, "09 ENA", "pacote ENA vazio (sem uViGs)")

    # ----------------------------- veredito -----------------------------
    c = rep.counts
    rep.raw("")
    rep.raw("-" * 64)
    rep.raw(f"  Resumo:  OK={c[OK]}  VAZIO-OK={c[EMPTY_OK]}  AVISO={c[WARN]}  FALHA={c[FAIL]}")
    if c[FAIL] > 0:
        verdict = ">>> VEREDITO: FALHA(S) DETECTADA(S) — investigar os itens [FALHA] acima."
    elif c[WARN] > 0:
        verdict = ">>> VEREDITO: INTEGRO COM AVISOS — etapas completaram; revise os [AVISO]."
    else:
        verdict = ">>> VEREDITO: PIPELINE INTEGRO — todas as etapas completaram e batem entre si."
    rep.raw("  " + verdict)
    rep.raw("-" * 64)

    text = "\n".join(rep.lines) + "\n"
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(text)
    print(text)

    if args.strict and c[FAIL] > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
