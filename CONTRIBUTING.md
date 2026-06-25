# Contributing

Thanks for your interest in improving this pipeline! Contributions of all sizes are
welcome — bug reports, documentation fixes, new options, or additional tools.

## Reporting bugs

Open a GitHub issue and include:

- The exact command you ran (`phages.py ...`).
- Operating system / environment (e.g., WSL2 Ubuntu 22.04, RAM, cores).
- The relevant log from `<output>/logs/` and the `RUN_HEALTHCHECK.txt` report.
- What you expected vs. what happened.

> Tip: the integrity check (`RUN_HEALTHCHECK.txt`, see README) usually pinpoints
> which stage failed and whether an empty result is legitimate or a real failure.

## Proposing changes

1. Fork the repository and create a branch (`git checkout -b my-fix`).
2. Keep the style of the surrounding code:
   - Snakemake rules in `workflow/rules/`, one stage per file.
   - Helper scripts in `workflow/scripts/` (standalone, argparse-based, no hidden state).
   - One isolated conda environment per tool in `workflow/envs/`.
3. Pin tool versions when a tool is known to break across releases (see the vRhyme
   pins in `workflow/envs/vrhyme.yaml` and `PIPELINE_DOCS.md` §10 for why).
4. Test a dry-run before submitting:
   ```bash
   python phages.py --mode genome --input <test>.fasta -o /tmp/test --block 1 -- -n
   ```
5. If you touched a stage, run the pipeline on a small sample and confirm the
   health check reports no `FALHA`.
6. Open a pull request describing the change and how you tested it.

## Adding a new tool / stage

- Add a `workflow/envs/<tool>.yaml` (channels `conda-forge`, `bioconda`).
- Add a `workflow/rules/<stage>.smk` with a clear log marker on completion
  (`echo "[<Tool>] OK -- ran to completion ..."`).
- Wire its expected outputs into `verify_run.py` so the integrity check covers it.

## Code of conduct

Be respectful and constructive. This is a research tool maintained by a small team;
please be patient with response times.
