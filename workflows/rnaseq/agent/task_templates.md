# Agent task templates

## Dry run

```bash
snakemake -np --use-conda
```

## Run workflow

```bash
snakemake --use-conda -j 32
```

## Validate sample sheet

```bash
python scripts/validate_samples.py --samples config/samples.tsv
```

Ensure `samples.tsv` includes the following columns:
- `sample`: Unique sample name.
- `fq1`: Path to R1 FASTQ or single-end FASTQ.
- `fq2`: Path to R2 FASTQ. Empty value means single-end sample.
Optional columns:
- `group`: Group or condition for differential expression analysis.
- `replicate`: Replicate ID for distinguishing biological replicates.

## Common debugging steps

1. Check `config/config.yml` paths.
2. Check `config/samples.tsv` and FASTQ existence.
3. Run `snakemake -np --use-conda`.
4. Inspect rule-specific logs under `results/logs/`.
5. Do not modify raw FASTQ files.
