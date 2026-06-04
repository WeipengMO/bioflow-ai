# lncRNA Discovery Agent Notes

## Before Running

- Read `README.md`, `agent/manifest.yml`, and `agent/io.contract.yml`.
- Deploy into a project directory with `scripts/deploy_pipeline.sh`; do not run production analyses directly inside `workflows/lncrna/`.
- Confirm `config/config.yml` points to a matched genome FASTA, annotation GTF, and STAR index from the same assembly.
- Validate `config/samples.tsv` with `python scripts/validate_samples.py --samples config/samples.tsv`.

## Standard Run

```bash
./run_snakemake.sh -n -p
./run_snakemake.sh
```

## Reporting Back

- Summarize candidate lncRNA files, expression matrices, BAM/bigWig outputs, and MultiQC report.
- If no candidates pass filtering, suggest checking `results/gffcompare/merged.stats` and relaxing `lncrna_filter` settings before changing upstream mapping parameters.
