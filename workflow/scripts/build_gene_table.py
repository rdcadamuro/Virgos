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


# PORTAO DE FAGO (taxonomia geNomad/ICTV). Os REALMS sao mistos (Duplodnaviria,
# Monodnaviria, Varidnaviria tem fagos E virus eucarioticos), por isso filtramos por
# CLASSE/FAMILIA. Lista comprehensiva de virus PROCARIOTICOS (bacteria + arqueia,
# DNA + RNA) -- nao so Caudoviricetes:
#   dsDNA caudados: Caudoviricetes (bact+arq) | ssDNA: Malgrandaviricetes (Microviridae),
#   Faserviricetes (Inoviridae) | ssRNA: Leviviricetes | dsRNA: Vidaverviricetes
#   (Cystoviridae) | arqueia: Tokiviricetes/Adnaviria + familias de arqueia.
DEFAULT_PHAGE_CLASSES = (
    "Caudoviricetes|Caudovirales|Leviviricetes|Faserviricetes|Vidaverviricetes|"
    "Malgrandaviricetes|Tokiviricetes|Adnaviria|"
    "Microviridae|Inoviridae|Plectroviridae|Tectiviridae|Corticoviridae|Autolykiviridae|"
    "Cystoviridae|Fiersviridae|Steitzviridae|Solspiviridae|Duinviridae|Plasmaviridae|"
    "Finnlakeviridae|Fuselloviridae|Halspiviridae|Thaspiviridae|Bicaudaviridae|"
    "Ampullaviridae|Clavaviridae|Globuloviridae|Guttaviridae|Spiraviridae|Ovaliviridae|"
    "Portogloboviridae|Turriviridae|Pleolipoviridae|Sphaerolipoviridae|Simuloviridae|"
    "Rudiviridae|Lipothrixviridae|Tristromaviridae"
)
# BLACKLIST EUCARIOTICA: classes/filos/familias que NUNCA sao fago (exclui mesmo se
# o virus ficar sub-classificado). Cobre retrovirus (HIV), herpes, virus gigantes,
# e os virus de RNA/ssDNA eucarioticos.
DEFAULT_EUK_CLASSES = (
    "Herviviricetes|Peploviricota|Nucleocytoviricota|Megaviricetes|Pokkesviricetes|"
    "Pimascovirales|Algavirales|Imitervirales|Pandoravirales|Naldaviricetes|"
    "Maveriviricetes|Adenoviridae|Lavidaviridae|Revtraviricetes|Artverviricota|"
    "Negarnaviricota|Cressdnaviricota|Pisoniviricetes|Stelpaviricetes|Alsuviricetes|"
    "Flasuviricetes|Magsaviricetes|Insthoviricetes|Monjiviricetes|Chunqiuviricetes|"
    "Ellioviricetes|Milneviricetes|Yunchangviricetes|Howeltoviricetes|Amabiliviricetes|"
    "Allassoviricetes|Miaviricetes|Repensiviricetes|Arfiviricetes|Mouviricetes|"
    "Quintoviricetes|Shotokuvirae"
)

_BIN_RE = re.compile(r"^(?P<sample>.+?)__(?:vamb|vrhyme)__bin\d+__(?P<contig>.+)$")


def parse_viral_regions(virus_summary, provirus):
    """Retorna:
      full_viral: {contig: taxonomy}      -- contigs totalmente virais (livres)
      prov: [(src, start, end, taxonomy)] -- proviruses (integrados = temperados)
    A taxonomia vem da coluna 'taxonomy' do virus_summary (inclui as linhas de provirus)."""
    full, prov = {}, []
    for r in read_tsv(virus_summary):
        sn = r.get("seq_name", "")
        tax = r.get("taxonomy", "") or ""
        if "|provirus_" in sn:
            m = re.search(r"\|provirus_(\d+)_(\d+)$", sn)
            src = sn.split("|provirus_")[0]
            if m:
                prov.append((src, int(m.group(1)), int(m.group(2)), tax))
        elif sn:
            full[sn] = tax
    # provirus.tsv reforça proviruses ausentes do virus_summary (taxonomia desconhecida)
    seen = {(s, a, b) for s, a, b, _ in prov}
    for r in read_tsv(provirus):
        src = r.get("source_seq") or r.get("seq_name", "")
        try:
            s, e = int(r.get("start", 0)), int(r.get("end", 0))
        except ValueError:
            continue
        src0 = src.split("|provirus_")[0]
        if src0 and e > 0 and (src0, s, e) not in seen:
            prov.append((src0, s, e, ""))
    return full, prov


def parse_phatyp(path, sample):
    """PhaTYP prediction -> {contig: TYPE} (temperate/virulent) SÓ desta amostra.
    A Accession do PhaTYP é o nome do vOTU (<sample>__<binner>__binN__<contig>);
    tiramos o prefixo p/ casar com o seq_name do geNomad (provirus -> contig de origem)."""
    out = {}
    for r in read_tsv(path):
        acc = r.get("Accession") or r.get("accession") or ""
        typ = (r.get("TYPE") or r.get("Type") or r.get("Pred") or "").strip().lower()
        if not acc or not typ:
            continue
        m = _BIN_RE.match(acc)
        if m:
            if m.group("sample") != sample:
                continue
            contig = m.group("contig").split("|provirus_")[0]
        else:
            contig = acc.split("|provirus_")[0]
        out[contig] = typ
    return out


def classify_phage(tax, phage_re, euk_re, lenient):
    """Decide se a taxonomia do geNomad indica FAGO (vírus procariótico). 3 vias:
      1) bate na whitelist de fago (classe/família procariótica) -> FAGO;
      2) bate na blacklist eucariótica (retro/herpes/gigantes/RNA euc.) -> NÃO-fago;
      3) sub-classificado (vírus sem rank útil): 'lenient' decide (default True =
         mantém, p/ NÃO perder fago novo/divergente; False = só whitelist).
    tax vazio (provirus reforçado pelo provirus.tsv) -> assume fago (prófago)."""
    if tax == "":
        return True
    if phage_re.search(tax):
        return True
    if euk_re.search(tax):
        return False
    return lenient


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
    ap.add_argument("--phatyp", default="", help="(opcional) PhaTYP prediction (lifestyle) -> filtro temperate")
    ap.add_argument("--phage-classes", default=DEFAULT_PHAGE_CLASSES,
                    help="whitelist (regex) de classes/famílias virais PROCARIÓTICAS (fagos)")
    ap.add_argument("--eukaryotic-classes", default=DEFAULT_EUK_CLASSES,
                    help="blacklist (regex) de táxons eucarióticos (nunca fago)")
    ap.add_argument("--phage-lenient", default="true",
                    help="true = mantém vírus sub-classificado (não perde fago novo); false = só whitelist")
    ap.add_argument("--viral-set", default="temperate_phage",
                    choices=["temperate_phage", "provirus", "all_viral"],
                    help="quais regiões virais entram como is_viral no enriquecimento")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mge_re = re.compile(args.mge_keywords, re.IGNORECASE)
    phage_re = re.compile(args.phage_classes, re.IGNORECASE)
    euk_re = re.compile(args.eukaryotic_classes, re.IGNORECASE)
    lenient = str(args.phage_lenient).strip().lower() in ("true", "1", "yes", "sim")
    full_viral, prov = parse_viral_regions(args.virus_summary, args.provirus)
    phatyp = parse_phatyp(args.phatyp, args.sample) if args.phatyp else {}
    mode = args.viral_set

    def is_phage(tax):
        return classify_phage(tax, phage_re, euk_re, lenient)

    # ----- seleciona o CONJUNTO VIRAL conforme o modo -----
    # provírus = integrados = TEMPERADOS (sempre); contigs livres precisam de
    # PhaTYP=temperate. Em ambos exige passar pelo PORTÃO DE FAGO (taxonomia).
    prov_sel, free_sel = [], set()
    n_prov_nonphage = n_free_nonphage = n_free_lytic = 0
    for (src, s, e, tax) in prov:
        if mode == "all_viral" or is_phage(tax):
            prov_sel.append((src, s, e))
        else:
            n_prov_nonphage += 1
    for contig, tax in full_viral.items():
        if mode == "all_viral":
            free_sel.add(contig)
        elif mode == "provirus":
            pass  # contigs livres não entram
        else:  # temperate_phage
            if not is_phage(tax):
                n_free_nonphage += 1
            elif phatyp.get(contig, "") != "temperate":
                n_free_lytic += 1
            else:
                free_sel.add(contig)
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

            viral = contig in free_sel or any(
                src == contig and not (ge < s or gs > e) for (src, s, e) in prov_sel)

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
    print(f"[build_gene_table] {args.sample}: {len(genes)} genes | viral_set={mode} "
          f"(regioes: {len(prov_sel)} provirus + {len(free_sel)} fago-livre-temperado) "
          f"-> {n_viral} genes virais")
    print(f"[build_gene_table] FILTRADOS do conjunto viral: "
          f"{n_prov_nonphage} provirus nao-fago, {n_free_nonphage} contigs nao-fago "
          f"(ex.: virus eucariotico/RNA), {n_free_lytic} fagos liticos (PhaTYP!=temperate)")
    print(f"[build_gene_table] {cat_str} -> {args.out}")


if __name__ == "__main__":
    main()
