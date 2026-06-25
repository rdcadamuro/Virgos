# Changelog

All notable changes to this project are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/); versions follow [SemVer](https://semver.org/).

## [0.1.0] — 2026-06-24

First public release. End-to-end validated on WSL2 (Ubuntu 22.04) with a complete
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

[0.1.0]: https://github.com/TODO-user/TODO-repo/releases/tag/v0.1.0
