# Virgos — Phage Identification & Analysis Pipeline

![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.0-brightgreen)
![License: MIT](https://img.shields.io/badge/license-MIT-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2-lightgrey)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20859599.svg)](https://doi.org/10.5281/zenodo.20859599)

> **Snakemake workflow** for recovery, QC, binning, annotation and abundance profiling of bacteriophages from complete genomes or metagenomes.
> Validated end-to-end on WSL2 Ubuntu 22.04 · 2026-06-24.

**Highlights**

- One launcher (`phages.py`), two modes: complete **genome** or **metagenome**.
- **One isolated conda environment per tool** — no dependency conflicts.
- **Whole-bin CheckV (N-linking)** — bin completeness as a single vMAG, per the vRhyme recommendation.
- **Automated integrity check** — catches *silent* failures (a tool that crashes but the run looks "successful").

---

## Overview

`phages` integrates six specialized tools into a single, reproducible pipeline with **isolated conda environments per tool** (no dependency conflicts). One command launcher (`phages.py`) handles both modes:

| Mode | Input | Steps |
|------|-------|-------|
| **`genome`** | Genome/assembly FASTA (± reads) | geNomad → CheckV → [binning if reads] → PhaBox2 → VITAP → summary |
| **`metagenome`** | Assembly + reads (Illumina/Nanopore) | geNomad → CheckV → binning → PhaBox2 → VITAP → **vOTU abundance** → summary |

---

## Workflow

```
FASTA [+ reads]
      │
      ▼
 01 geNomad          viral contig identification (language model + MMseqs2)
      │
      ▼
 02 CheckV           QC: completeness, contamination, viral genes
      │
      ▼
 03 vRhyme + VAMB    co-coverage binning → filter → ANI/AF dereplication → vOTUs
      │                                                         │
      ▼                                                         ▼ (metagenome only)
 04 PhaBox2          lifestyle (PhaTYP), taxonomy (PhaMer+PhaGCN),
    CHERRY host      host prediction via CRISPR/alignment
      │
      ▼
 05 VITAP            ICTV-framework taxonomy via bipartite graph
      │
      ▼
 06 Summary          viral_bins_table.tsv + per_sample_summary.tsv
      │
      ▼
 ✓ Integrity check   RUN_HEALTHCHECK.txt (scans logs + cross-checks stages)
```

> When reads are present, each bin's scaffolds are joined with an N-spacer and
> CheckV is run on the **whole bin** (N-linking) for an accurate vMAG completeness
> estimate. See [Quality control & run integrity](#quality-control--run-integrity).

---

## Tools

| Tool | Version | Purpose |
|------|---------|---------|
| [geNomad](https://github.com/apcamargo/genomad) | ≥ 1.7 | Viral sequence identification |
| [CheckV](https://bitbucket.org/berkeleylab/checkv) | ≥ 1.0.1 | Genome quality assessment |
| [vRhyme](https://github.com/AnantharamanLab/vRhyme) | v1.1.0 (pinned env, see note) | Co-coverage viral binning |
| [VAMB](https://github.com/RasmussenLab/vamb) | ≥ 4.0 | VAE-based metagenomic binning |
| [PhaBox2](https://github.com/KennthShang/PhaBOX) | ≥ 2.0 | Lifestyle, taxonomy, host |
| [VITAP](https://github.com/DrKaiyangZheng/VITAP) | = 1.10 | ICTV taxonomic assignment |
| [CoverM](https://github.com/wwood/CoverM) | ≥ 0.7 | vOTU abundance profiling |
| [Snakemake](https://snakemake.github.io/) | ≥ 9.0 | Workflow orchestration |

Each tool runs in its own isolated conda environment — no version conflicts.

---

## Installation

**Requirements:** Linux (Ubuntu 22.04 tested), conda/miniforge, ≥ 8 GB RAM, ≥ 20 GB disk.

> **WSL2 note:** output directory must be on a Linux filesystem (`/home/...`), not `/mnt/c/` (NTFS breaks MMseqs2 memory mapping).

```bash
# 1. Install Miniforge (if not already installed)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b && ~/miniforge3/bin/conda init bash && exec bash

# 2. Clone the repository
git clone https://github.com/<username>/<repo>.git ~/viral_pipeline
cd ~/viral_pipeline

# 3. Create the orchestration environment (from the provided file)
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate viralpipe

# 4. Download databases and build all per-tool conda environments (Block 0)
python phages.py --block 0
```

After Block 0, edit `config/config.yaml` to set database paths, then run.

> See [PIPELINE_DOCS.md](PIPELINE_DOCS.md) for full installation details, database setup, and known compatibility fixes.

---

## Quick Start

```bash
conda activate viralpipe
cd ~/viral_pipeline

# Complete genome + paired-end reads
python phages.py \
  --mode genome \
  --input /path/to/genome.fasta \
  --reads1 /path/to/R1.fastq.gz \
  --reads2 /path/to/R2.fastq.gz \
  --read-type ilu-paired \
  -o /path/to/output \
  --block 1

# Metagenome — multiple samples via TSV table
python phages.py \
  --mode metagenome \
  --samples config/samples.tsv \
  -o /path/to/output \
  --block 1

# Interactive mode (prompts for all options)
python phages.py
```

### Sample table (`samples.tsv`) — metagenome mode

```tsv
sample	assembly	reads1	reads2	read_type
SAMPLE_A	/data/A.fasta	/data/A_R1.fq.gz	/data/A_R2.fq.gz	ilu-paired
SAMPLE_B	/data/B.fasta	/data/B_R1.fq.gz		ilu-single
SAMPLE_C	/data/C.fasta	/data/C_nano.fq.gz		nano
```

`read_type` options: `ilu-paired`, `ilu-single`, `nano`.

---

## Outputs

Results are written to numbered subdirectories inside `<output>/`:

```
<output>/
├── 01_genomad/           viral contigs (geNomad)
├── 02_checkv/            per-contig CheckV; <sample>_bins/ = whole-bin CheckV (N-linking)
├── 03_bins_vrhyme_vamb/  bins, linked_bins.fna, dereplicated vOTUs (votus.fna)
├── 04_phabox2/           lifestyle, taxonomy, host prediction
├── 05_vitap/             ICTV taxonomic lineage
├── 06_votu_abundance/    [metagenome] TPM matrix (vOTU × sample)
├── 06_summary_final/ (genome) / 07_summary_final/ (metagenome)
│   ├── viral_bins_table.tsv     ← main result: one row per vOTU
│   ├── per_sample_summary.tsv   ← counts by quality/lifestyle per sample
│   └── RUN_HEALTHCHECK.txt       ← integrity report (see below)
└── logs/                  per-rule logs, each with a completion marker
```

### Main output columns (`viral_bins_table.tsv`)

| Column | Description |
|--------|-------------|
| `votu_rep` | vOTU representative contig ID |
| `rep_sample` | sample of origin |
| `binner` | `vamb` or `vrhyme` |
| `checkv_quality` | Complete / High-quality / Medium-quality / Low-quality |
| `completeness` | CheckV estimated completeness (%) |
| `lifestyle` | `lysogenic`, `lytic`, or `-` (PhaTYP) |
| `taxonomy` | VITAP lineage (empty if below confidence threshold) |
| `mean_tpm` | [metagenome] mean abundance across samples in TPM |
| `prevalence` | [metagenome] number of samples where detected |

---

## Quality control & run integrity

Two features make results trustworthy and reduce the risk of **silent failures**
(a tool that crashes mid-run while the pipeline still reports "done"):

**1. Whole-bin CheckV (N-linking).** When reads are present, a bin (vMAG) may contain
several scaffolds that are fragments of the same genome. Following the vRhyme
recommendation, the pipeline joins a bin's scaffolds with a 1000-N spacer and runs
CheckV on the **whole bin**, so completeness reflects the vMAG rather than the
longest fragment. (Sequences sent downstream to dereplication/annotation are the
real scaffolds — the N-spacer is used only for the CheckV estimate.)

**2. Automated integrity check.** Every full run ends with `verify_run.py`
(rule `healthcheck`), which writes `RUN_HEALTHCHECK.txt` and prints it to the
console. It detects problems two independent ways:

- **Log scan** — fatal signatures (`Traceback`, `No module named`, `command not
  found`, `core dumped`, out-of-memory, etc.).
- **Cross-stage consistency** — e.g. "geNomad found N>0 contigs but CheckV evaluated
  0 rows" → failure; "N vOTUs but 0 PhaBox2 predictions" → failure.

Each item is labelled `OK` / `EMPTY-OK` (legitimately empty) / `WARN` / `FAIL`,
with a final verdict. Run it standalone anytime:

```bash
python workflow/scripts/verify_run.py \
    --outdir <output> --mode genome --samples <output>/samples.tsv \
    --out <output>/06_summary_final/RUN_HEALTHCHECK.txt
```

Add `--strict` to make a `FAIL` exit non-zero (aborts the pipeline).

---

## Functional annotation: resistome, mobilome & enrichment

An optional block (`0X_annotation/`) characterises the **resistome** (ARGs) and
**mobilome** (mobile genetic elements) of the recovered vOTUs and tests, by
**permutation**, whether they are *enriched* in the viral fraction vs the whole genome.

- **Resistome** — ABRicate across **all bundled databases** (`abricate_dbs: "all"`
  → NCBI AMR, CARD, ResFinder, ARG-ANNOT, MEGARES, VFDB, Victors, EcOH, ecoli_vf,
  upec_expec_vf, PlasmidFinder, BacMet2), run with `--threads`, each DB's raw hits
  kept for provenance. Hits are deduplicated **by locus** (reciprocal-overlap +
  union-find: the same gene reported by several DBs collapses to one locus listing
  the agreeing DBs, without chaining adjacent genes). Outputs a **per-DB summary**
  (`resistome_by_db.tsv`) and a self-verifying **dedup checkpoint**
  (`resistome_dedup_checkpoint.json`, read by the integrity check). The bundled DBs
  need no download; block 0 records a version **manifest** for reproducibility.
- **Mobilome** — consolidates several tools, each emitting coordinate calls:
  geNomad (provirus / conjugation / integrase), **IntegronFinder** (integrons),
  **MobileElementFinder** (IS/transposons), **CONJScan/MacSyFinder** (conjugative
  systems — the installable, host-agnostic alternative to ICEfinder), **tRNAscan-SE**.
  ICEfinder is supported as an optional manual hook.
- **2-layer deduplication** (`mge_dedup.py`) — layer 1 collapses redundancy *within*
  each tool; layer 2 merges overlapping calls *across* tools into non-redundant
  regions, tagged with which tools agree (agreement = confidence). The **reciprocal-overlap**
  criterion merges the same element reported by different tools but keeps *nested*
  elements distinct (an IS inside an ICE stays separate). A self-verifying
  **checkpoint** asserts no residual overlaps (read by the integrity check).
- **Per-category enrichment** — `build_gene_table.py` labels every genome gene by
  category (`is_arg`, `is_conjugation`, `is_integrase`, `is_transposase`, `is_integron`,
  `is_conjugative_system`, `is_IS`, `is_mge`); `permutation_enrichment.py` runs a
  **stratified permutation** (+ hypergeometric cross-check) per category, with
  guardrails (`NO_BACKGROUND` for pure-phage input, `UNDERPOWERED` for rare events).
  Each category is labelled by the **same tool applied to the whole genome** so the
  test is valid: geNomad (arg/conjugation/integrase/transposase/recombinase),
  **IntegronFinder** (`is_integron`), **MobileElementFinder** (`is_IS`) and
  **CONJScan** (`is_conjugative_system`) all run a second time on the full assembly —
  not only geNomad. Without this, `is_IS`/`is_conjugative_system` would have no
  non-viral background and stay permanently `UNDERPOWERED`.
- **Reports** — `annotation_master.tsv` (per-vOTU detail) + `annotation_stats.log`
  (statistical summary with automatic ENRICHED/DEPLETED/n.s./UNDERPOWERED interpretation).

> The enrichment test only applies when the input genome has a **non-viral fraction**
> (a bacterial genome / MAG with prophages). For a complete phage genome it reports
> `NO_BACKGROUND`. See [PIPELINE_DOCS.md §13](PIPELINE_DOCS.md#13-anotação-funcional-resistoma--mobiloma--enriquecimento).

---

## ENA submission package (MIUViG)

A modular, mode-aware block (`ena.smk`) auto-builds a ready-to-submit **ENA**
package for the recovered uViGs, following the **MIUViG** standard (Roux et al.
2019) and ENA's **MAG assembly** route (Webin-CLI). It runs at the end
(`stage_ena`, included in the full run when `ena.enable: true`) and writes to
`0X_summary_final/ena_submission/`:

- `ena_uvig_metadata.tsv` — master table with **every MIUViG field per uViG**,
  auto-filled from what the pipeline already produces: `assembly_qual` /
  completeness / contamination (**CheckV**), `source_uvig` & `detec_type` (provirus
  vs free virus), taxonomy (**VITAP**), OTU clustering (ANI/AF), gene prediction
  & viral identification (**geNomad**), `COVERAGE` (CoverM or contig `_cov_`),
  `PLATFORM` (from read type).
- `ena_sample_registration.tsv` + `ena_samples.xml` — one **SAMPLE per uViG**
  (each MAG needs its own sample) with the MIUViG attributes and checklist tag.
- `manifests/<alias>.manifest` — a **Webin-CLI** MAG manifest per uViG.
- `fasta/<alias>.fasta.gz` — each uViG sequence with a clean header.
- `submit_all.sh` — loops `webin-cli` (validate → submit) over all manifests.
- `ENA_SUBMISSION_GUIDE.md` — step-by-step guide citing all software + versions.

Only the fields that **only the submitter knows** are left as `TODO_PREENCHER`
(collection date, location, nucleic-acid extraction, upstream assembler, study/
sample accessions, tax_id) — and these can be pre-set in `config.yaml → ena:`.
The integrity check reports an `09 ENA` item (uViGs = manifests; TODO count).
See [PIPELINE_DOCS.md §14](PIPELINE_DOCS.md#14-submissão-ao-ena-pacote-miuvig-dos-uvigs).

---

## Configuration

Key parameters in [`config/config.yaml`](config/config.yaml):

```yaml
databases:
  genomad: "/home/<user>/viral_dbs/genomad_db"
  checkv:  "/home/<user>/viral_dbs/checkv-db-v1.5"
  phabox:  "/home/<user>/viral_dbs/phabox_db_v2"
  vitap:   "/home/<user>/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR"

genomad:
  min_virus_score: 0.7
  extra: "--cleanup --splits 4"   # --splits 4 required on WSL2 with 8 GB RAM

derep:
  min_ani: 95                     # vOTU clustering: ANI >= 95%
  min_af:  85                     # vOTU clustering: alignment fraction >= 85%
```

---

## Known Issues & Compatibility Notes

| Issue | Fix |
|-------|-----|
| **vRhyme silently crashed mid-run** (`pkg_resources` removed in setuptools≥81; pickled models need scikit-learn<1.3) | Env pins `python<3.12`, `scikit-learn<1.3`, `setuptools<81`; error masking removed (PIPELINE_DOCS §10.8) |
| geNomad OOM on WSL2 (8 GB RAM) | Add `--splits 4` to `genomad.extra` in config.yaml |
| MMseqs2 fails on NTFS (`/mnt/c/`) | Set output to Linux filesystem (`/home/...`) |
| PhaBox2 >= 2.0 + pandas >= 2.0 incompatibility | Durable, auto-applied source patch (`patch_phabox2.py`) |
| VITAP Diamond DB version mismatch (Figshare) | Rebuild `.dmnd` with installed Diamond version (`diamond makedb`) |
| VITAP `uniref90.dmnd` not included in DB package | Pipeline auto-creates empty fallback file to skip UniRef90 BLAST step |
| VITAP `db_dir` at wrong level after ZIP extraction | Point config to `ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/` (the subdirectory) |

Full details and patch code: [PIPELINE_DOCS.md](PIPELINE_DOCS.md#10-problemas-conhecidos-e-correções-aplicadas).

---

## Validation

Tested end-to-end on sample **MS13934** (complete Caudoviricetes phage genome + paired-end Illumina reads):

| Step | Status | Notes |
|------|--------|-------|
| geNomad | ✔ | 5 viral contigs identified |
| CheckV | ✔ | 1 Medium-quality (57.3%), 4 Low-quality |
| vRhyme + VAMB | ✔ | vRhyme **ran to completion** (0 bins, legitimate — single sample) + 5 VAMB bins → 5 vOTUs |
| PhaBox2 | ✔ | bin1: lysogenic, *Enterococcus termitis* host (CRISPR); bin2: lytic; all Caudoviricetes |
| VITAP | ✔ | Ran to completion; no lineage below Caudoviricetes (expected for <60% complete fragments) |
| Summary | ✔ | viral_bins_table.tsv + per_sample_summary.tsv generated |
| Integrity check | ✔ | `RUN_HEALTHCHECK.txt` → **INTEGRO COM AVISOS** (0 failures) |

Platform: WSL2 Ubuntu 22.04 · 8 GB RAM · 8 cores · Snakemake 9.23.1

> The integrity check was also verified to **catch a planted silent failure**
> (`FAIL` verdict), so a clean report is meaningful — not a rubber stamp.

---

## Related pipelines

This workflow follows the now-standard viral-metagenomics pattern (viral identification
→ QC → vOTU clustering → taxonomy → abundance). Established pipelines in the same space:

| Pipeline | Engine | Notes |
|----------|--------|-------|
| [VirMake](https://github.com/Rounge-lab/VirMake) | Snakemake | QC→assembly→identification→vOTU→taxonomy/function→abundance |
| [MVP](https://github.com/RasmussenLab/mvp) (Modular Viromics Pipeline) | — | identify, filter, cluster, annotate, bin |
| [Hecatomb](https://github.com/shandley/hecatomb) | Snakemake | integrated viral metagenomics platform |
| [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) | Nextflow/CWL | EBI MGnify viral detection |
| [niaid/virome-pipeline](https://github.com/niaid/virome-pipeline) | Snakemake | geNomad/CheckV → vOTUs |

**Where this pipeline fits:** it is intentionally focused and reproducible rather than
all-encompassing — a single launcher, isolated per-tool conda environments, a dual
genome/metagenome mode, the vRhyme-style **whole-bin CheckV**, and an **automated
integrity check** that hardens it against the silent tool failures that are common
when chaining many bioinformatics tools. See [`awesome-virome`](https://github.com/shandley/awesome-virome)
for a broader catalogue.

---

## Citation

**To cite this pipeline**, use the metadata in [`CITATION.cff`](CITATION.cff) (GitHub
shows a *"Cite this repository"* button) or the Zenodo DOI:

- **All versions (concept DOI, recommended):** [10.5281/zenodo.20859599](https://doi.org/10.5281/zenodo.20859599)
- **This release (v1.0.0):** [10.5281/zenodo.20859600](https://doi.org/10.5281/zenodo.20859600)

If you use this pipeline in your research, **please also cite the underlying tools**:

- **geNomad:** Camargo et al., *Nature Biotechnology* 2023. https://doi.org/10.1038/s41587-023-01953-y
- **CheckV:** Nayfach et al., *Nature Biotechnology* 2021. https://doi.org/10.1038/s41587-020-00774-7
- **vRhyme:** Kieft et al., *Nucleic Acids Research* 2022. https://doi.org/10.1093/nar/gkac341
- **VAMB:** Nissen et al., *Nature Biotechnology* 2021. https://doi.org/10.1038/s41587-020-00777-4
- **PhaBox2:** Shang et al., "PhaBOX2: an enhanced web server for discovering and analyzing viral contigs in metagenomic data", *Nucleic Acids Research* 2026 (advance access). https://doi.org/10.1093/nar/gkag382
- **VITAP:** Zheng et al., "VITAP: a high precision tool for DNA and RNA viral classification based on meta-omic data", *Nature Communications* 16, 2226 (2025). https://doi.org/10.1038/s41467-025-57500-7
- **Diamond:** Buchfink et al., *Nature Methods* 2021. https://doi.org/10.1038/s41592-021-01101-x
- **Snakemake:** Mölder et al., *F1000Research* 2021. https://doi.org/10.12688/f1000research.29032.2

---

## Documentation & project files

| File | What it is |
|------|-----------|
| [PIPELINE_DOCS.md](PIPELINE_DOCS.md) | Full technical reference: architecture, every output, all known issues & fixes, the vRhyme bug chain, N-linking, the integrity check |
| [CHANGELOG.md](CHANGELOG.md) | Version history |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to report bugs and contribute |
| [CITATION.cff](CITATION.cff) | Citation metadata (fill in author/DOI before publishing) |
| [config/config.yaml](config/config.yaml) | All tunable parameters and database paths |

> Note: `PIPELINE_DOCS.md` is currently written in Portuguese; an English translation
> is planned. This README is the canonical English entry point.

---

## License

MIT — see [LICENSE](LICENSE).

