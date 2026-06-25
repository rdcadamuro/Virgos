# Changelog

All notable changes to this project are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/); versions follow [SemVer](https://semver.org/).

## [1.0.0] — 2026-06-25

First stable release (**Virgos**). Adds a full functional-analysis layer and an
ENA submission module on top of the validated core. End-to-end validated on WSL2
(Ubuntu 22.04) with a complete bacterial genome carrying a prophage (sample
MS13934, *Enterococcus*); integrity check reports `INTEGRO` (0 failures).

### Added
- **Functional annotation block** (`annotation.smk`, `enrichment.smk`):
  - **Resistome:** ABRicate across **all bundled databases** (`abricate_dbs: "all"`
    → NCBI, CARD, ResFinder, ARG-ANNOT, MEGARES, VFDB, Victors, EcOH, ecoli_vf,
    upec_expec_vf, PlasmidFinder, BacMet2), with `--threads`, per-DB raw outputs,
    a per-DB summary, and a block-0 version manifest.
  - **Mobilome:** consolidates geNomad + IntegronFinder + MobileElementFinder +
    CONJScan + tRNAscan per vOTU.
  - **2-layer deduplication** (`mge_dedup.py`, `dedup_abricate.py`): reciprocal-
    overlap + union-find, with self-verifying **checkpoints** (residual overlaps
    must be 0) read by the integrity check.
  - **Per-category permutation enrichment** of ARG/MGE in the viral fraction vs the
    whole genome (stratified permutation + hypergeometric cross-check + guardrails),
    with whole-genome IntegronFinder / MobileElementFinder / CONJScan feeding the
    `is_integron` / `is_IS` / `is_conjugative_system` categories.
- **ENA submission package** (`ena.smk`, `build_ena_tables.py`): auto-builds the
  **MIUViG** metadata table, per-uViG Webin-CLI **MAG manifests**, sample XML +
  registration sheet, gzipped FASTAs, a batch `submit_all.sh`, and a step-by-step
  guide — auto-filled from CheckV (quality/completeness/contamination), VITAP
  (taxonomy), derep (ANI/AF), geNomad, coverage and platform. Mode-aware (genome
  or metagenome).
- Integrity check extended with `06 resistoma`/`06 dedup-ARG`, per-category
  `06 enriquecimento`, `06 dedup-MGE`, and `09 ENA` items.

### Fixed
- **CONJScan was non-functional** everywhere: `--models CONJScan` is invalid in
  CONJScan 2.x (needs `CONJScan/Chromosome` + `CONJScan/Plasmids`) and the install
  guard checked `macsydata available` (remote list) instead of installed models →
  it always degraded to empty. Fixed model paths + `--user` install + dir check.
- **MobileElementFinder `pkg_resources` crash** (setuptools≥81) — pinned `setuptools<81`.
- **Enrichment only used geNomad:** `is_IS` / `is_conjugative_system` were structurally
  empty; now MobileElementFinder and CONJScan run on the whole genome and feed them.
- Corrected tool citations (**PhaBox2** DOI pointed to an unrelated paper; **VITAP** DOI missing).

## [0.1.0] — 2026-06-24

Core pipeline. End-to-end validated on WSL2 (Ubuntu 22.04) with a complete
phage genome (sample MS13934, *Caudoviricetes*).

### Added
- Dual-mode launcher `phages.py`: `genome` and `metagenome`.
- Pipeline stages: geNomad → CheckV → vRhyme+VAMB binning → dereplication (vOTUs)
  → PhaBox2 (lifestyle/taxonomy/host) → VITAP (ICTV taxonomy) → (metagenome) CoverM
  abundance → integrated summary tables.
- Isolated conda environment per tool (`--use-conda`), no dependency conflicts.
- **Whole-bin CheckV via N-linking** (vRhyme recommendation): scaffolds of each bin
  are joined with an N-spacer before CheckV so completeness reflects the whole vMAG.
- **Intelligent integrity check** (`verify_run.py` + `healthcheck` rule): scans logs
  for fatal signatures and cross-checks stage-to-stage consistency to catch silent
  failures; writes `RUN_HEALTHCHECK.txt`.
- Per-stage completion markers in every log ("ran to completion" vs "0 results").

### Fixed
- **vRhyme silent failure (critical):** the env resolved to `setuptools>=81`
  (no `pkg_resources`) and `scikit-learn>=1.3` (incompatible with vRhyme's pickled
  models), so vRhyme crashed mid-run while the rule masked it with `|| true`. Pinned
  `python<3.12`, `scikit-learn<1.3`, `setuptools<81`; removed the error masking and
  added real error-vs-zero-bins detection.
- **PhaBox2 + pandas≥2 incompatibility:** durable, idempotent source patch
  (`patch_phabox2.py`) re-applied automatically by the `phabox` rule.
- **VITAP:** Diamond DB rebuilt for the installed Diamond version; auto-fallback when
  `uniref90.dmnd` is absent; `db_dir` pointed to the correct extracted subfolder.
- **geNomad on WSL2:** `--splits 4` to avoid MMseqs2 OOM; output must be on a Linux
  filesystem (NTFS breaks MMseqs2 memory mapping).

[1.0.0]: https://github.com/rdcadamuro/Virgos/releases/tag/v1.0.0
[0.1.0]: https://github.com/rdcadamuro/Virgos/releases/tag/v0.1.0
