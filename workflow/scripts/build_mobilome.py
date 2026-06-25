#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tabela de MOBILOMA por vOTU, a partir das anotacoes que o geNomad JA produz
(sem rodar ferramenta nova).

Para cada vOTU (representante dereplicado), consolida:
  - topology         : DTR / ITR / Provirus / No terminal repeats (virus_summary)
  - is_provirus      : se o contig e um provirus integrado (geNomad find_proviruses)
  - n_amr_genes      : genes com anotacao AMR do geNomad (annotation_amr)
  - amr_genes        : descricoes desses genes
  - n_conjugation    : genes com marcador de conjugacao (annotation_conjscan)
  - n_integrase / n_transposase / n_other_mge : por palavra-chave na descricao

Mapeamento vOTU -> contig do geNomad: tira o prefixo de bin (<sample>__<binner>__
bin<N>__) quando presente; o resto e o seq_name do geNomad (inclui sufixo
|provirus_<s>_<e> dos provirus).
"""
import argparse
import csv
import os
import re

NA = {"", "NA", "na", "-", "nan", "None"}


def read_fasta_names(path):
    names = []
    if not os.path.exists(path):
        return names
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0].strip())
    return names


def read_fasta_lengths(path):
    """-> {name: length} (name = 1o token do header)."""
    lengths, name, n = {}, None, 0
    if not os.path.exists(path):
        return lengths
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0].strip(); n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


# offset p/ traduzir coord do frame geNomad -> frame do vOTU.
# Para um provirus extraido (...|provirus_<s>_<e>), o geNomad numera os genes no
# frame do contig de origem (comeca em <s>); a sequencia do vOTU comeca em 1.
# offset = <s> - 1. Para vOTU = contig inteiro, offset = 0.
_PROV_SUFFIX = re.compile(r"\|provirus_(\d+)_(\d+)$")


def votu_offset(contig):
    m = _PROV_SUFFIX.search(contig)
    return (int(m.group(1)) - 1) if m else 0


def read_tsv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path, errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]
    return header, rows


def gene_to_seq(gene_id):
    """gene id do geNomad = <seq_name>_<idx> -> seq_name (idx e o ultimo _N)."""
    return gene_id.rsplit("_", 1)[0] if "_" in gene_id else gene_id


def load_sample_index(genomad_dir, sample, mge_re):
    base = os.path.join(genomad_dir, sample)
    summ = os.path.join(base, f"{sample}_summary")
    genes_fp = os.path.join(summ, f"{sample}_virus_genes.tsv")
    vsum_fp = os.path.join(summ, f"{sample}_virus_summary.tsv")
    prov_fp = os.path.join(base, f"{sample}_find_proviruses", f"{sample}_provirus.tsv")

    # topology por seq_name
    topo = {}
    _, vrows = read_tsv(vsum_fp)
    for r in vrows:
        topo[r.get("seq_name", "")] = r.get("topology", "")

    # genes por seq_name -> contadores + COORDENADAS dos genes MGE (p/ element calls).
    # mge_coords[seq] = [(start, end, subtype), ...] no FRAME do seq_name do geNomad
    # (p/ provirus, esse frame e o do contig de origem: ex. 131..14917).
    genes = {}
    mge_coords = {}
    _, grows = read_tsv(genes_fp)
    for r in grows:
        seq = gene_to_seq(r.get("gene", ""))
        d = genes.setdefault(seq, {"amr": [], "conj": 0, "integ": 0, "transp": 0, "mge": 0})
        amr = r.get("annotation_amr", "")
        if amr not in NA:
            d["amr"].append(amr)
        is_conj = r.get("annotation_conjscan", "") not in NA
        if is_conj:
            d["conj"] += 1
        desc = (r.get("annotation_description", "") or "")
        is_integ = is_transp = False
        if mge_re.search(desc):
            d["mge"] += 1
            low = desc.lower()
            if "integrase" in low:
                d["integ"] += 1; is_integ = True
            elif "transposase" in low:
                d["transp"] += 1; is_transp = True
        # registra coordenada do gene se for sinal de mobiloma com posicao util
        subtype = ("integrase" if is_integ else "conjugation" if is_conj
                   else "transposase" if is_transp else None)
        if subtype is not None:
            try:
                gs, ge = int(r.get("start", 0)), int(r.get("end", 0))
            except (ValueError, TypeError):
                gs = ge = 0
            if ge > 0:
                mge_coords.setdefault(seq, []).append((min(gs, ge), max(gs, ge), subtype))

    # provirus: seq_names + span (start/end do contig de origem) -> element call
    prov = set()
    prov_span = {}  # seq_name -> (src_contig, start, end)
    _, prows = read_tsv(prov_fp)
    for r in prows:
        sn = r.get("seq_name", "")
        if sn:
            prov.add(sn)
            try:
                s, e = int(r.get("start", 0)), int(r.get("end", 0))
            except (ValueError, TypeError):
                s = e = 0
            prov_span[sn] = (r.get("source_seq", ""), s, e)

    return {"topo": topo, "genes": genes, "prov": prov,
            "mge_coords": mge_coords, "prov_span": prov_span}


BIN_RE = re.compile(r"^(?P<sample>.+?)__(?:vamb|vrhyme)__bin\d+__(?P<contig>.+)$")


def parse_votu(rep, samples):
    """-> (sample_or_None, contig_geNomad)."""
    m = BIN_RE.match(rep)
    if m:
        return m.group("sample"), m.group("contig")
    # sem prefixo de bin (modo genoma sem reads): rep == contig do geNomad
    return (samples[0] if len(samples) == 1 else None), rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--votus", required=True)
    ap.add_argument("--genomad-dir", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--mge-keywords",
                    default="integrase|transposase|recombinase|relaxase|resolvase|mobiliz|conjugat")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-calls", default="",
                    help="(opcional) TSV de element calls normalizados no frame do vOTU")
    args = ap.parse_args()

    mge_re = re.compile(args.mge_keywords, re.IGNORECASE)
    idx = {s: load_sample_index(args.genomad_dir, s, mge_re) for s in args.samples}

    reps = read_fasta_names(args.votus)
    lengths = read_fasta_lengths(args.votus)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    calls = []  # element calls no frame do vOTU (seq=rep)

    def add_call(rep, start, end, etype, score):
        L = lengths.get(rep)
        s = max(1, int(start)); e = int(end)
        if L:
            e = min(e, L)
        if e < s:
            s, e = e, s
        if e >= 1:
            calls.append((rep, s, e, etype, score))

    cols = ["votu_rep", "sample", "topology", "is_provirus", "n_amr_genes",
            "amr_genes", "n_conjugation_genes", "n_integrase", "n_transposase",
            "n_mge_total"]
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for rep in reps:
            samp, contig = parse_votu(rep, args.samples)
            cands = [samp] if samp in idx else list(idx)
            chosen, g, topo, isprov, ix_hit = "", None, "", False, None
            for s in cands:
                ix = idx[s]
                if contig in ix["genes"] or contig in ix["topo"] or contig in ix["prov"]:
                    chosen = s
                    g = ix["genes"].get(contig)
                    topo = ix["topo"].get(contig, "")
                    isprov = (contig in ix["prov"]) or ("provirus" in contig.lower()) \
                        or (topo.lower() == "provirus")
                    ix_hit = ix
                    break
            if g is None:
                g = {"amr": [], "conj": 0, "integ": 0, "transp": 0, "mge": 0}
                isprov = isprov or ("provirus" in contig.lower())
            out.write("\t".join(str(x) for x in [
                rep, chosen, topo, "Yes" if isprov else "No",
                len(g["amr"]), ";".join(g["amr"]) if g["amr"] else "-",
                g["conj"], g["integ"], g["transp"], g["mge"],
            ]) + "\n")

            # ---- element calls (frame do vOTU) ----
            if ix_hit is None:
                continue
            off = votu_offset(contig)  # geNomad coord -> vOTU coord = coord - off
            # 1) provirus: a sequencia inteira do vOTU e o elemento (1..len)
            if isprov:
                L = lengths.get(rep)
                if L:
                    add_call(rep, 1, L, "provirus", 1)
                else:
                    span = ix_hit["prov_span"].get(contig)
                    if span and span[2] > 0:
                        add_call(rep, span[1] - off, span[2] - off, "provirus", 1)
            # 2) genes MGE (integrase/conjugacao/transposase) com coordenada
            for (gs, ge, subtype) in ix_hit["mge_coords"].get(contig, []):
                etype = {"integrase": "integrase", "conjugation": "conjugative_system",
                         "transposase": "transposase"}.get(subtype, "mge_gene")
                add_call(rep, gs - off, ge - off, etype, 1)

    if args.out_calls:
        os.makedirs(os.path.dirname(args.out_calls) or ".", exist_ok=True)
        with open(args.out_calls, "w", newline="") as oc:
            w = csv.writer(oc, delimiter="\t")
            w.writerow(["seq", "start", "end", "element_type", "tool", "score"])
            for (rep, s, e, etype, score) in sorted(calls, key=lambda x: (x[0], x[1])):
                w.writerow([rep, s, e, etype, "geNomad", f"{float(score):.2f}"])
        print(f"[build_mobilome] element calls (geNomad): {len(calls)} -> {args.out_calls}")

    print(f"[build_mobilome] {len(reps)} vOTUs -> {args.out}")


if __name__ == "__main__":
    main()
