#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gera o PACOTE DE SUBMISSAO ENA dos uViGs (Uncultivated Virus Genomes).

Auto-alimenta TUDO que o pipeline ja deriva, seguindo o padrao MIUViG (Minimum
Information about an Uncultivated Virus Genome, Roux et al. 2019, GSC) e a rota
de submissao de assembly do ENA (Webin-CLI, tipo MAG). Campos que so o submissor
sabe (data de coleta, local, extracao de acido nucleico) ficam como TODO claro.

Entradas (todas opcionais -> degrada com nota):
  --per-bin      viral_bins_table.tsv  (espinha: 1 linha por vOTU; quality,
                 completeness, lifestyle, taxonomia VITAP, binner, amostras)
  --checkv-bins  1+ quality_summary.tsv do CheckV nos BINS (miuvig_quality,
                 contamination, completeness_method) — chave: contig_id=rep_unit
  --vitap        vitap_taxonomy.tsv (lineage + Confidence_level)
  --mobilome     mobilome_table.tsv (is_provirus, topology)
  --samples      samples.tsv (read_type -> PLATFORM)
  --votus        votus.fna (sequencias -> FASTA por vOTU; cov no nome)
  --abundance    (metagenoma) votu_abundance_tpm.tsv p/ COVERAGE
  --config-json  secao 'ena' do config + versoes/parametros (geNomad, CheckV...)

Saidas em --outdir:
  ena_uvig_metadata.tsv        master legivel: TODOS os campos MIUViG (auto+TODO)
  ena_sample_registration.tsv  1 linha por uViG p/ registrar SAMPLE no Webin
  ena_samples.xml              SAMPLE_SET XML pronto (checklist do config)
  manifests/<alias>.manifest   manifesto Webin-CLI (ASSEMBLY tipo MAG) por uViG
  fasta/<alias>.fasta.gz       sequencia do uViG (header limpo)
  submit_all.sh                helper: loop chamando webin-cli em cada manifesto
  ENA_SUBMISSION_GUIDE.md      passo-a-passo + softwares/versoes citados
"""
import argparse
import csv
import gzip
import json
import os
import re

NA = {"", "NA", "na", "-", "nan", "None", None}
TODO = "TODO_PREENCHER"


# ----------------------------- leitura ------------------------------
def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, errors="ignore", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    """-> {name: seq} (name = 1o token do header)."""
    seqs, name, chunks = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
    with open(path, errors="ignore") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = ln[1:].split()[0].strip()
                chunks = []
            else:
                chunks.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


# ----------------------------- mapeamentos --------------------------
_COV_RE = re.compile(r"_cov_(\d+\.?\d*)", re.I)
_PROV_RE = re.compile(r"\|provirus_\d+_\d+$")


def parse_cov_from_name(name):
    m = _COV_RE.search(name)
    return f"{float(m.group(1)):.1f}" if m else ""


def platform_from_readtype(rt):
    rt = (rt or "").lower()
    if rt.startswith("ilu"):
        return "ILLUMINA"
    if rt == "nano":
        return "OXFORD_NANOPORE"
    return ""


# CheckV miuvig_quality -> vocabulario controlado MIUViG assembly_qual
def assembly_qual(miuvig_q, checkv_q):
    m = (miuvig_q or "").lower()
    c = (checkv_q or "").lower()
    if "high-quality" in m or "complete" in c:
        if "complete" in c:
            return "Finished genome"
        return "High-quality draft genome (single contig; >=90% complete; <5% contamination)"
    return "Genome fragment(s)"


def topology_to_struc(topology):
    """pred_genome_struc do MIUViG = segmented/non-segmented (default non-segmented
    p/ a maioria dos fagos). Mantemos a topologia (DTR/ITR/provirus) numa coluna a parte."""
    return "non-segmented"


def sanitize_alias(sample, idx):
    s = re.sub(r"[^A-Za-z0-9]+", "_", sample or "S").strip("_")
    return f"{s}_vOTU{idx:03d}"


# ----------------------------- core ---------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-bin", required=True)
    ap.add_argument("--checkv-bins", nargs="*", default=[])
    ap.add_argument("--vitap", default="")
    ap.add_argument("--mobilome", default="")
    ap.add_argument("--samples", default="")
    ap.add_argument("--votus", default="")
    ap.add_argument("--abundance", default="")
    ap.add_argument("--config-json", default="")
    ap.add_argument("--mode", default="metagenome")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    cfg = {}
    if args.config_json and os.path.exists(args.config_json):
        try:
            cfg = json.load(open(args.config_json))
        except (ValueError, OSError):
            cfg = {}

    def c(key, default=""):
        v = cfg.get(key, default)
        return v if v not in NA else default

    os.makedirs(args.outdir, exist_ok=True)
    fasta_dir = os.path.join(args.outdir, "fasta")
    manif_dir = os.path.join(args.outdir, "manifests")
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(manif_dir, exist_ok=True)

    per_bin = read_tsv(args.per_bin)

    # CheckV bins (varios arquivos) -> {contig_id: row}
    checkv = {}
    for f in args.checkv_bins:
        for r in read_tsv(f):
            cid = r.get("contig_id", "")
            if cid:
                checkv[cid] = r

    # VITAP -> {genome_id: (lineage, confidence)}
    vitap = {}
    for r in read_tsv(args.vitap):
        gid = r.get("Genome_ID") or r.get("genome_id") or ""
        if gid:
            vitap[gid] = (r.get("lineage", ""), r.get("Confidence_level", ""))

    # mobiloma -> {votu_rep: (is_provirus, topology)}
    mob = {}
    for r in read_tsv(args.mobilome):
        rep = r.get("votu_rep", "")
        if rep:
            mob[rep] = (r.get("is_provirus", "No"), r.get("topology", ""))

    # samples -> {sample: platform}
    plat_by_sample = {}
    for r in read_tsv(args.samples):
        plat_by_sample[r.get("sample", "")] = platform_from_readtype(r.get("read_type", ""))

    # abundancia (metagenoma) -> {votu_rep: mean_coverage} se a tabela tiver
    abund = {}
    for r in read_tsv(args.abundance):
        key = r.get("votu_rep") or r.get("votu") or r.get("contig") or ""
        cov = r.get("mean_coverage") or r.get("coverage") or r.get("Mean") or ""
        if key and cov not in NA:
            abund[key] = cov

    seqs = read_fasta(args.votus)

    # ---- defaults do config (com fallback p/ TODO onde so humano sabe) ----
    genomad_ver = c("genomad_version", "1.12.0")
    checkv_ver = c("checkv_version", "1.1.1")
    min_score = c("genomad_min_virus_score", "0.7")
    ani = c("derep_ani", "95")
    af = c("derep_af", "85")
    source_default = c("source_uvig", "viral fraction metagenome (virome)")
    assembly_sw = c("assembly_software", "") or f"{TODO} (assembler; versao; parametros)"
    plat_override = c("sequencing_platform", "")
    molecule = c("molecule_type", "genomic DNA")
    enrich = c("virus_enrich_appr", "") or TODO
    nucl_ext = c("nucl_acid_ext", "") or TODO
    coll_date = c("collection_date", "") or TODO
    geo_loc = c("geographic_location", "") or TODO
    pred_type = c("pred_genome_type", "uncharacterized")
    study = c("study_accession", "") or "TODO_STUDY (PRJEBxxxxx)"
    center = c("center_name", "")
    checklist = c("sample_checklist", "ERC000011")
    mag_taxid = c("uvig_tax_id", "")  # taxid p/ os uViGs (ex.: 12429 unclassified viruses)
    mag_taxname = c("uvig_scientific_name", "uncultured virus")

    feat_pred = f"geNomad/pyrodigal-gv; {genomad_ver}; default"
    vir_ident = f"geNomad; {genomad_ver}; end-to-end, min_virus_score={min_score}"
    ref_db = f"geNomad-DB; CheckV-DB (v{checkv_ver} marker set)"
    sim_search = f"geNomad marker search (MMseqs2 profiles); {genomad_ver}"
    compl_appr_base = f"CheckV; {checkv_ver}"
    otu_appr = f"{ani}% ANI; {af}% AF; greedy incremental (CheckV anicalc/aniclust)"
    tax_class_method = f"VITAP; {c('vitap_version','latest')}; best-hit LCA vs ICTV/RefSeq/IMG-VR"

    MIUVIG_COLS = [
        "votu_id", "ena_alias", "sample_source", "sequence_length", "number_contig",
        "source_uvig", "detec_type", "vir_ident_software",
        "assembly_qual", "checkv_quality", "completeness", "compl_score", "compl_appr",
        "contamination", "pred_genome_type", "pred_genome_struc", "topology",
        "feat_pred", "ref_db", "sim_search_meth", "tax_class",
        "vitap_lineage", "vitap_confidence", "otu_class_appr", "lifestyle",
        "host_pred_appr", "sort_tech", "virus_enrich_appr", "nucl_acid_ext",
        "collection_date", "geographic_location", "coverage", "sequencing_platform",
        "assembly_software", "molecule_type",
    ]

    rows = []
    samples_rows = []
    n_provirus = 0
    for idx, b in enumerate(per_bin, 1):
        rep = b.get("votu_rep", "")
        unit = b.get("rep_unit", "")
        samp = b.get("rep_sample", "") or (b.get("all_samples", "") or "").split(";")[0]
        all_s = b.get("all_samples", "") or samp
        alias = sanitize_alias(samp, idx)
        cv = checkv.get(unit, {})
        is_prov = (mob.get(rep, ("No", ""))[0] == "Yes") or (cv.get("provirus", "No") == "Yes") \
            or bool(_PROV_RE.search(rep))
        topology = mob.get(rep, ("", ""))[1] or ("Provirus" if is_prov else "")
        if is_prov:
            n_provirus += 1
        source_uvig = "provirus (UViG)" if is_prov else source_default
        detec = "provirus (UViG)" if is_prov else "independent sequence (UViG)"

        miuvig_q = cv.get("miuvig_quality", "")
        checkv_q = b.get("checkv_quality", "") or cv.get("checkv_quality", "")
        compl = b.get("completeness", "") or cv.get("completeness", "")
        contam = cv.get("contamination", "0")
        compl_method = cv.get("completeness_method", "") or "AAI-based"
        seq = seqs.get(rep, "")
        seqlen = len(seq) if seq else (cv.get("contig_length", "") or "")

        # taxonomia: per_bin.taxonomy (ja vem do VITAP) + vitap direto p/ confianca
        tax_lineage = b.get("taxonomy", "") or vitap.get(rep, ("", ""))[0] \
            or vitap.get(unit, ("", ""))[0]
        tax_conf = vitap.get(rep, ("", ""))[1] or vitap.get(unit, ("", ""))[1]
        tax_lineage = tax_lineage if tax_lineage not in NA else "unclassified (sem hit VITAP)"

        # coverage: metagenoma->abundancia; senao parse do nome (_cov_)
        cov = ""
        if args.mode == "metagenome":
            cov = abund.get(rep, "") or abund.get(unit, "")
        if cov in NA:
            cov = parse_cov_from_name(rep)
        if cov in NA:
            cov = TODO

        # plataforma: override do config, senao do read_type das amostras
        plats = sorted({plat_by_sample.get(s, "") for s in all_s.split(";") if s})
        plats = [p for p in plats if p]
        platform = plat_override or (",".join(plats) if plats else TODO)

        compl_str = f"{compl}%" if compl not in NA else ""
        compl_score = (f"high;{compl_str}" if _to_float(compl) >= 90 else
                       f"med;{compl_str}" if _to_float(compl) >= 50 else
                       f"low;{compl_str}" if compl not in NA else "")

        row = {
            "votu_id": rep, "ena_alias": alias, "sample_source": all_s,
            "sequence_length": seqlen, "number_contig": "1",
            "source_uvig": source_uvig, "detec_type": detec,
            "vir_ident_software": vir_ident,
            "assembly_qual": assembly_qual(miuvig_q, checkv_q),
            "checkv_quality": checkv_q, "completeness": compl_str,
            "compl_score": compl_score, "compl_appr": f"{compl_appr_base}; {compl_method}",
            "contamination": f"{contam}%" if contam not in NA else "",
            "pred_genome_type": pred_type, "pred_genome_struc": topology_to_struc(topology),
            "topology": topology or "No terminal repeats",
            "feat_pred": feat_pred, "ref_db": ref_db, "sim_search_meth": sim_search,
            "tax_class": tax_class_method, "vitap_lineage": tax_lineage,
            "vitap_confidence": tax_conf or "NA", "otu_class_appr": otu_appr,
            "lifestyle": b.get("lifestyle", ""), "host_pred_appr": "not provided",
            "sort_tech": "not applicable", "virus_enrich_appr": enrich,
            "nucl_acid_ext": nucl_ext, "collection_date": coll_date,
            "geographic_location": geo_loc, "coverage": cov,
            "sequencing_platform": platform, "assembly_software": assembly_sw,
            "molecule_type": molecule,
        }
        rows.append(row)

        # registro de SAMPLE (1 por uViG / MAG)
        samples_rows.append({
            "alias": alias, "tax_id": mag_taxid or TODO,
            "scientific_name": mag_taxname,
            "sample_title": f"uViG {alias} ({source_uvig})",
            "sample_description": f"Uncultivated virus genome derivado de {samp} "
                                  f"(binner={b.get('binner','')}); {row['assembly_qual']}; "
                                  f"completeness={compl_str}; taxonomia(VITAP)={tax_lineage}",
            "source_uvig": source_uvig, "assembly_qual": row["assembly_qual"],
            "completeness_score": compl_score, "number_contig": "1",
            "collection_date": coll_date, "geographic_location": geo_loc,
            "isolation_source": source_uvig, "tax_class": tax_lineage,
        })

        # FASTA por uViG (header = alias limpo) + manifesto
        if seq:
            with gzip.open(os.path.join(fasta_dir, f"{alias}.fasta.gz"), "wt") as fa:
                fa.write(f">{alias}\n")
                for i in range(0, len(seq), 70):
                    fa.write(seq[i:i + 70] + "\n")
        write_manifest(manif_dir, alias, study, molecule, assembly_sw, platform,
                       cov, b.get("rep_sample", "") or samp, rep)

    # ---------------- escreve as tabelas ----------------
    _write_tsv(os.path.join(args.outdir, "ena_uvig_metadata.tsv"), MIUVIG_COLS, rows)
    sreg_cols = ["alias", "tax_id", "scientific_name", "sample_title", "sample_description",
                 "source_uvig", "assembly_qual", "completeness_score", "number_contig",
                 "collection_date", "geographic_location", "isolation_source", "tax_class"]
    _write_tsv(os.path.join(args.outdir, "ena_sample_registration.tsv"), sreg_cols, samples_rows)
    write_sample_xml(os.path.join(args.outdir, "ena_samples.xml"), samples_rows,
                     checklist, center, mag_taxid, mag_taxname)
    write_submit_helper(os.path.join(args.outdir, "submit_all.sh"), [r["ena_alias"] for r in rows])
    write_guide(os.path.join(args.outdir, "ENA_SUBMISSION_GUIDE.md"), rows, n_provirus,
                checklist, study, genomad_ver, checkv_ver, ani, af, args.mode)

    n_todo = sum(1 for r in rows for k in ("collection_date", "geographic_location",
                 "virus_enrich_appr", "nucl_acid_ext") if str(r[k]).startswith("TODO"))
    print(f"[build_ena_tables] {len(rows)} uViGs -> {args.outdir} "
          f"(provirus={n_provirus}; campos TODO a preencher por uViG: "
          f"{n_todo // max(1, len(rows))} de 4)")


def _to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return -1.0


def _write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_manifest(manif_dir, alias, study, molecule, program, platform, coverage,
                   run_sample, votu_id):
    """Manifesto Webin-CLI (context=genome, ASSEMBLY_TYPE=MAG)."""
    cov = coverage if coverage and not str(coverage).startswith("TODO") else "1"
    lines = [
        f"STUDY\t{study}",
        f"SAMPLE\tTODO_SAMPLE_ACC (apos registrar {alias})",
        f"ASSEMBLYNAME\t{alias}",
        "ASSEMBLY_TYPE\tMetagenome-Assembled Genome (MAG)",
        f"COVERAGE\t{cov}",
        f"PROGRAM\t{program}",
        f"PLATFORM\t{platform}",
        f"MOLECULETYPE\t{molecule}",
        f"DESCRIPTION\tuViG {alias} (origem: {votu_id})",
        f"FASTA\tfasta/{alias}.fasta.gz",
    ]
    with open(os.path.join(manif_dir, f"{alias}.manifest"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _xml_escape(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def write_sample_xml(path, samples_rows, checklist, center, taxid, taxname):
    """SAMPLE_SET XML pronto p/ submeter via Webin (REST/programatico)."""
    out = ['<?xml version="1.0" encoding="UTF-8"?>', "<SAMPLE_SET>"]
    attr_keys = ["source_uvig", "assembly_qual", "completeness_score", "number_contig",
                 "collection_date", "geographic_location", "isolation_source"]
    for s in samples_rows:
        center_attr = f' center_name="{_xml_escape(center)}"' if center else ""
        out.append(f'  <SAMPLE alias="{_xml_escape(s["alias"])}"{center_attr}>')
        out.append(f'    <TITLE>{_xml_escape(s["sample_title"])}</TITLE>')
        out.append("    <SAMPLE_NAME>")
        out.append(f"      <TAXON_ID>{_xml_escape(s.get('tax_id') or taxid or 'TODO')}</TAXON_ID>")
        out.append(f"      <SCIENTIFIC_NAME>{_xml_escape(taxname)}</SCIENTIFIC_NAME>")
        out.append("    </SAMPLE_NAME>")
        out.append(f"    <DESCRIPTION>{_xml_escape(s['sample_description'])}</DESCRIPTION>")
        out.append("    <SAMPLE_ATTRIBUTES>")
        # checklist obrigatorio (ENA usa o atributo ENA-CHECKLIST)
        out.append("      <SAMPLE_ATTRIBUTE>")
        out.append("        <TAG>ENA-CHECKLIST</TAG>")
        out.append(f"        <VALUE>{_xml_escape(checklist)}</VALUE>")
        out.append("      </SAMPLE_ATTRIBUTE>")
        for k in attr_keys:
            out.append("      <SAMPLE_ATTRIBUTE>")
            out.append(f"        <TAG>{k}</TAG>")
            out.append(f"        <VALUE>{_xml_escape(s.get(k, 'TODO'))}</VALUE>")
            out.append("      </SAMPLE_ATTRIBUTE>")
        out.append("    </SAMPLE_ATTRIBUTES>")
        out.append("  </SAMPLE>")
    out.append("</SAMPLE_SET>")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def write_submit_helper(path, aliases):
    """Helper bash: loop chamando webin-cli em cada manifesto (validate -> submit)."""
    lines = [
        "#!/bin/bash",
        "# Submissao em lote dos uViGs ao ENA via Webin-CLI.",
        "# Pre-requisitos: estudo + samples JA registrados; preencher STUDY/SAMPLE nos .manifest.",
        "#   1) export WEBIN_USER=Webin-XXXXX ; export WEBIN_PASS=...",
        "#   2) baixe webin-cli.jar: https://github.com/enasequence/webin-cli/releases",
        "set -uo pipefail",
        'JAR="${WEBIN_CLI_JAR:-webin-cli.jar}"',
        'MODE="${1:-validate}"   # validate | submit',
        "for m in manifests/*.manifest; do",
        '  echo "=== $m ($MODE) ==="',
        '  java -jar "$JAR" -context=genome -manifest="$m" \\',
        '       -userName="$WEBIN_USER" -password="$WEBIN_PASS" -"$MODE" || echo "  FALHOU: $m"',
        "done",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_guide(path, rows, n_provirus, checklist, study, gver, cver, ani, af, mode):
    n = len(rows)
    auto = ("source_uvig, detec_type, vir_ident_software, assembly_qual, checkv_quality, "
            "completeness, compl_score, compl_appr, contamination, pred_genome_struc, "
            "feat_pred, ref_db, sim_search_meth, tax_class, vitap_lineage, otu_class_appr, "
            "number_contig, sequence_length, coverage, sequencing_platform, molecule_type")
    todo = ("collection_date, geographic_location, virus_enrich_appr, nucl_acid_ext, "
            "assembly_software (assembler a montante), study_accession, sample_accession, tax_id")
    g = f"""# Guia de submissao ao ENA — uViGs (Uncultivated Virus Genomes)

Gerado automaticamente pelo pipeline. **{n} uViGs** ({n_provirus} proviruses),
modo `{mode}`. Padrao **MIUViG** (Roux et al. 2019, Nat Biotechnol) + rota de
**assembly MAG** do ENA (Webin-CLI).

## Arquivos deste pacote
| Arquivo | Para que serve |
|---|---|
| `ena_uvig_metadata.tsv` | Tabela-mestre legivel: TODOS os campos MIUViG por uViG (auto + TODO). |
| `ena_sample_registration.tsv` | 1 linha por uViG p/ registrar o SAMPLE (planilha Webin). |
| `ena_samples.xml` | SAMPLE_SET pronto p/ submissao programatica (checklist `{checklist}`). |
| `manifests/<alias>.manifest` | Manifesto Webin-CLI (ASSEMBLY_TYPE = MAG) por uViG. |
| `fasta/<alias>.fasta.gz` | Sequencia do uViG (header limpo) p/ o Webin-CLI. |
| `submit_all.sh` | Loop que chama o webin-cli em cada manifesto (validate/submit). |

## Campos preenchidos automaticamente
{auto}.

## Campos que VOCE precisa preencher (so o submissor sabe)
{todo}.
Procure `TODO_PREENCHER` em `ena_uvig_metadata.tsv` / `*.manifest` / `ena_samples.xml`.

## Passo a passo
1. **Conta Webin**: crie/login em https://www.ebi.ac.uk/ena/submit/webin/ (`Webin-XXXXX`).
2. **Estudo (Project)**: registre 1 estudo -> obtem `PRJEBxxxxx`. Substitua `{study}`
   em todos os `manifests/*.manifest` (campo STUDY).
3. **Checklist do sample**: confirme o checklist de UViG/ambiental adequado no
   navegador de checklists do ENA (https://www.ebi.ac.uk/ena/browser/checklists).
   O default usado aqui foi `{checklist}` — ajuste se o ENA indicar outro (ex.: GSC MIxS).
4. **tax_id**: uViGs sem classificacao -> use um taxid de virus nao classificado
   (ex.: *unclassified viruses*, taxid 12429) ou o sugerido pela linhagem VITAP em
   `vitap_lineage`. Preencha `uvig_tax_id` no config p/ propagar automaticamente.
5. **Registrar samples**: suba `ena_samples.xml` (1 SAMPLE por uViG) via Webin REST
   ou cole `ena_sample_registration.tsv` na planilha interativa do Webin. Guarde os
   accessions `ERSxxxxxxx` retornados e preencha o campo SAMPLE nos `*.manifest`.
6. **Submeter assemblies**: `bash submit_all.sh validate` (checa) e depois
   `bash submit_all.sh submit`. Cada uViG vira 1 assembly MAG.

## Softwares citados (proveniencia MIUViG)
- Identificacao viral: **geNomad {gver}** (end-to-end).
- Qualidade/completeness/contaminacao: **CheckV {cver}** (assembly_qual MIUViG = `miuvig_quality`).
- Binning: **vRhyme + VAMB**; vOTUs por **{ani}% ANI / {af}% AF** (CheckV anicalc/aniclust).
- Taxonomia: **VITAP** (coluna `vitap_lineage`).
- Predicao de genes/anotacao: **geNomad / pyrodigal-gv**.

> Observacao: o ENA exige que cada MAG tenha um SAMPLE proprio derivado da amostra
> de origem ("sample derived from"). Os samples gerados aqui ja sao 1-por-uViG.
"""
    with open(path, "w") as fh:
        fh.write(g)


if __name__ == "__main__":
    main()
