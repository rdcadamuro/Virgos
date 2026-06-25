#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONJScan/MacSyFinder por vOTU -> sistemas conjugativos (T4SS/relaxase) com COORDENADA.

Alternativa instalavel e host-agnostica ao ICEfinder (CONJScan e o modulo adotado
pelo ICEfinder 2.0). Para cada contig do votus.fna:
  1) prodigal -> proteinas em ORDEM no genoma (.faa) + coordenadas (.gff/sco);
  2) macsyfinder --models CONJScan all --db-type ordered_replicon --sequence-db <faa>;
  3) parse do best_solution.tsv -> genes que compoem cada sistema;
  4) SPAN do sistema = (min start, max end) das coordenadas dos genes do sistema,
     no FRAME do contig do vOTU (prodigal numera os genes na ordem -> usamos o GFF).

Aplica CAMADA 1 de dedup (multi-modelo no MESMO locus) via mge_dedup.reciprocal_overlap
e emite ELEMENT CALLS normalizados:
    seq, start, end, element_type=conjugative_system, tool=CONJScan, score

Degradacao graciosa: sem macsyfinder/prodigal, sem CONJScan instalado, ou sem
vOTUs -> escreve TSV vazio COM HEADER + nota; NUNCA falha dura (exit 0 sempre).
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

# importa o motor de dedup ja validado (mesmo diretorio scripts/)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from mge_dedup import reciprocal_overlap
except Exception:  # fallback defensivo: dedup vira no-op se o import falhar
    def reciprocal_overlap(a, b):
        return 0.0

CALL_HDR = ["seq", "start", "end", "element_type", "tool", "score"]


def log(msg, fh):
    fh.write(msg + "\n")
    fh.flush()


def read_fasta(path):
    """Itera (name, seq) do FASTA (name = 1o token do header)."""
    name, chunks = None, []
    with open(path, errors="ignore") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = ln[1:].split()[0].strip()
                chunks = []
            else:
                chunks.append(ln.strip())
    if name is not None:
        yield name, "".join(chunks)


def run_prodigal(seq, contig_name, workdir, logfh):
    """Roda prodigal num contig isolado. Retorna (faa_path, {prot_id: (start,end)})
    ou (None, {}) em falha. prot_id segue a numeracao do prodigal: <contig>_<idx>."""
    fna = os.path.join(workdir, "c.fna")
    faa = os.path.join(workdir, "c.faa")
    gff = os.path.join(workdir, "c.gff")
    with open(fna, "w") as fh:
        fh.write(f">{contig_name}\n{seq}\n")
    # -p meta: robusto p/ contigs curtos/virais; -f gff p/ coordenadas
    try:
        subprocess.run(["prodigal", "-i", fna, "-a", faa, "-f", "gff", "-o", gff,
                        "-p", "meta", "-q"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        log(f"  [prodigal] falhou em {contig_name}: {e}", logfh)
        return None, {}
    coords = {}
    idx = 0
    with open(gff, errors="ignore") as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] != "CDS":
                continue
            idx += 1
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            # prodigal escreve o id da proteina como ">_<idx>" no faa header;
            # casamos pela ORDEM (idx 1..N), que e como o macsyfinder ve em ordered_replicon
            coords[idx] = (min(s, e), max(s, e))
    return (faa if os.path.exists(faa) and os.path.getsize(faa) > 0 else None), coords


def parse_macsy_best(out_dir):
    """Le best_solution.tsv (ordered) ou all_systems.tsv (fallback) -> lista de
    (system_id, [hit_id,...]). hit_id = id da proteina do faa (prodigal:
    <contig>_<idx>). Trata cabecalhos repetidos (linha que comeca com 'replicon')."""
    cand = [os.path.join(out_dir, "best_solution.tsv"),
            os.path.join(out_dir, "all_best_solutions.tsv"),
            os.path.join(out_dir, "all_systems.tsv")]
    path = next((p for p in cand if os.path.exists(p) and os.path.getsize(p) > 0), None)
    if not path:
        return []
    systems = {}
    idx = None
    with open(path, errors="ignore") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"):
                continue
            f = ln.split("\t")
            if f and f[0] == "replicon":          # cabecalho (pode repetir)
                idx = {c: i for i, c in enumerate(f)}
                continue
            if idx is None:
                continue

            def gv(col):
                i = idx.get(col)
                return f[i] if i is not None and i < len(f) else ""
            sys_id = gv("sys_id") or gv("model_fqn") or "sys"
            hit = gv("hit_id")
            if hit:
                systems.setdefault(sys_id, []).append(hit)
    return list(systems.items())


def hit_to_idx(hit_id):
    """prodigal protein id <contig>_<idx> -> idx (int) ou None."""
    try:
        return int(hit_id.rsplit("_", 1)[-1])
    except (ValueError, IndexError):
        return None


def dedup_layer1(calls, overlap):
    """Funde calls do MESMO contig cujos spans se sobrepoem reciprocamente >= overlap
    (camada 1: multi-modelo CONJScan no mesmo locus). O(n^2), n pequeno por contig."""
    by_seq = {}
    for c in calls:
        by_seq.setdefault(c["seq"], []).append(c)
    out = []
    for seq, cs in by_seq.items():
        used = [False] * len(cs)
        for i in range(len(cs)):
            if used[i]:
                continue
            s, e, sc = cs[i]["start"], cs[i]["end"], cs[i]["score"]
            used[i] = True
            for j in range(i + 1, len(cs)):
                if used[j]:
                    continue
                if reciprocal_overlap((s, e), (cs[j]["start"], cs[j]["end"])) >= overlap:
                    s = min(s, cs[j]["start"]); e = max(e, cs[j]["end"])
                    sc = max(sc, cs[j]["score"]); used[j] = True
            out.append({"seq": seq, "start": s, "end": e,
                        "element_type": "conjugative_system", "tool": "CONJScan",
                        "score": sc})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--votus", required=True)
    ap.add_argument("--out", required=True, help="TSV de element calls (CONJScan)")
    ap.add_argument("--models", default="CONJScan/Chromosome CONJScan/Plasmids",
                    help="conjuntos de modelos (separados por espaco); CONJScan 2.x "
                         "tem subconjuntos Chromosome e Plasmids")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--log", default="")
    args = ap.parse_args()

    logfh = open(args.log, "a") if args.log else sys.stderr
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    def finish(calls):
        with open(args.out, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(CALL_HDR)
            for c in sorted(calls, key=lambda x: (x["seq"], x["start"])):
                w.writerow([c["seq"], c["start"], c["end"], c["element_type"],
                            c["tool"], f"{float(c['score']):.2f}"])
        log(f"[CONJScan] element calls escritos: {len(calls)} -> {args.out}", logfh)

    if not os.path.exists(args.votus) or os.path.getsize(args.votus) == 0:
        log("[CONJScan] sem vOTUs -> vazio. NAO e erro.", logfh)
        finish([]); return 0
    if not shutil.which("macsyfinder") or not shutil.which("prodigal"):
        log("[CONJScan] macsyfinder/prodigal ausentes -> vazio. NAO e erro.", logfh)
        finish([]); return 0

    base_work = args.workdir or tempfile.mkdtemp(prefix="conjscan_")
    os.makedirs(base_work, exist_ok=True)

    raw_calls = []
    n_contigs = 0
    for name, seq in read_fasta(args.votus):
        if not seq:
            continue
        n_contigs += 1
        cwork = os.path.join(base_work, f"c{n_contigs}")
        if os.path.isdir(cwork):
            shutil.rmtree(cwork, ignore_errors=True)
        os.makedirs(cwork, exist_ok=True)
        faa, coords = run_prodigal(seq, name, cwork, logfh)
        if not faa or not coords:
            log(f"  [CONJScan] {name}: sem proteinas (prodigal) -> pulado", logfh)
            continue
        # 1 macsyfinder por conjunto de modelos (Chromosome + Plasmids)
        for ms in args.models.split():
            out_dir = os.path.join(cwork, "macsy_" + ms.replace("/", "_"))
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)
            cmd = ["macsyfinder", "--models", ms, "all",
                   "--db-type", "ordered_replicon", "--sequence-db", faa,
                   "-o", out_dir, "-w", "1"]
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                log(f"  [CONJScan] {name}/{ms}: macsyfinder retornou erro ({e}); seguindo", logfh)
                continue
            for sys_id, hits in parse_macsy_best(out_dir):
                spans = [coords[i] for i in (hit_to_idx(h) for h in hits)
                         if i is not None and i in coords]
                if not spans:
                    continue
                s = min(a for a, _ in spans)
                e = max(b for _, b in spans)
                raw_calls.append({"seq": name, "start": s, "end": e, "score": len(spans)})
                log(f"  [CONJScan] {name}: {ms} sistema {sys_id} -> {s}-{e} ({len(spans)} genes)", logfh)

    # CAMADA 1: multi-modelo no mesmo locus -> 1 call
    deduped = dedup_layer1(raw_calls, args.overlap)
    log(f"[CONJScan] contigs={n_contigs} raw_systems={len(raw_calls)} "
        f"-> calls(camada1)={len(deduped)}", logfh)
    finish(deduped)
    if not args.workdir:
        shutil.rmtree(base_work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
