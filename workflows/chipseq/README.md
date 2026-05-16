# ChIP-seq / CUT&Tag Workflow

## Overview

This BioFlowAI workflow processes paired-end or single-end ChIP-seq/CUT&Tag FASTQ files into cleaned FASTQs, aligned duplicate-removed BAMs, normalized bigWig tracks, deepTools profiles, MACS3 peaks, replicate consensus peaks, ChIPQC reports, and HOMER motif enrichment results.

The workflow is config-driven and self-contained. It avoids project-specific hard-coded paths in rules and keeps reusable workflow logic under `rules/`, environments under `envs/`, helper scripts under `scripts/`, and agent metadata under `agent/`.

## Inputs

Required files:

```text
config/config.yml
config/samples.yml
raw_data/*.fastq.gz or raw_data/*.fq.gz
```

Optional files:

```text
config/replicates.yml
config/fastq_manifest.tsv
config/chipqc_sample_sheet.tsv
```

For paired-end data, automatic FASTQ discovery supports common names such as:

```text
sample_R1.fastq.gz / sample_R2.fastq.gz
sample_1.fq.gz     / sample_2.fq.gz
sample.R1.fq.gz    / sample.R2.fq.gz
sample-R1.fq.gz    / sample-R2.fq.gz
sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
```

For irregular file names or multi-lane FASTQs that were not merged, provide `config/fastq_manifest.tsv` and set `fastq_manifest` in `config/config.yml`.

## Configuration

Start from the examples:

```bash
cp config/config.example.yml config/config.yml
cp config/samples.example.yml config/samples.yml
cp config/replicates.example.yml config/replicates.yml
```

Important config keys:

| Key | Meaning |
|---|---|
| `mode` | `pe` or `se` |
| `raw_data_dir` | Directory containing FASTQ files |
| `fastq_manifest` | Optional TSV with `sample`, `read1`, `read2` columns |
| `genome` | Bowtie2 index prefix passed to `bowtie2 -x` |
| `regions_bed` | BED regions for deepTools `computeMatrix` |
| `gsize` | MACS3 genome size such as `hs`, `mm`, or a numeric effective genome size |
| `samples` | YAML file defining treatment and control samples |
| `replicates` | YAML file defining replicate groups |
| `call_peak_modes` | `with_control`, `without_control`, or both |
| `call_peak_types` | `narrow`, `broad`, or both |
| `enable_chipqc` | Generate ChIPQC report |
| `enable_homer` | Run HOMER motif enrichment |

HOMER parallelism is controlled through `threads.homer` in `config/config.yml`.

## Samples and controls

`config/samples.yml` separates treatment samples from control strategy.

Supported control strategies:

- `none`: call peaks without control only.
- `matched`: each treatment has a matched control in `controls.matched`.
- `pooled`: merge `controls.pooled_samples` into `controls.pooled_name` and use it for all treatments.

## Outputs

Main outputs:

```text
aligned_data/{sample}.sorted.rmdup.bam
aligned_data/{sample}.sorted.rmdup.bam.bai
tracks/{sample}.sorted.rmdup.CPM.bw
deeptools_profile/{sample}.scale.png
macs3_results/{narrow,broad,narrow_no_control,broad_no_control}/
replicate_intersect/{group}_intersect.bed
reports/chipqc/
reports/homer/
```

Logs are written to:

```text
logs/
```

## Quick Start

Recommended: deploy a clean project directory first.

```bash
cd workflows/chipseq
./scripts/deploy_pipeline.sh /path/to/chipseq_project
cd /path/to/chipseq_project
# put FASTQ files under raw_data/
./run_snakemake.sh -n --cores 1
```

If `/path/to/chipseq_project` already contains files and you want to overwrite deploy-managed files and links:

```bash
cd workflows/chipseq
./scripts/deploy_pipeline.sh /path/to/chipseq_project --force
```

Manual setup is still supported:

```bash
cd workflows/chipseq
cp config/config.example.yml config/config.yml
cp config/samples.example.yml config/samples.yml
cp config/replicates.example.yml config/replicates.yml
mkdir -p raw_data
# put FASTQ files under raw_data/
snakemake -np --use-conda
snakemake --use-conda -j 16
```

A convenience wrapper is also provided:

```bash
./run_snakemake.sh -n
./run_snakemake.sh
```

## FASTQ pair detection

The previous workflow inferred read 2 by slicing `sep`, for example `_1` to `_2`. This is fragile for mixed naming schemes and requires rebuilding paths from sample names.

The refactored workflow scans `raw_data_dir` once at parse time and builds an in-memory sample-to-FASTQ index. Rules then retrieve read paths from that index. This is faster, easier to validate, and supports several common naming schemes. For non-standard names, use an explicit manifest:

```tsv
sample	read1	read2
WT_rep1	raw_data/a_L001_R1.fastq.gz	raw_data/a_L001_R2.fastq.gz
```

You can generate a draft manifest:

```bash
python scripts/make_fastq_manifest.py --raw-data-dir raw_data --output config/fastq_manifest.tsv
```

Then set:

```yaml
fastq_manifest: config/fastq_manifest.tsv
```

## Troubleshooting

Run a dry-run first:

```bash
snakemake -np --use-conda
```

Common checks:

- FASTQ sample names match `config/samples.yml`.
- `genome` is a Bowtie2 index prefix, not just a FASTA path.
- `regions_bed` matches the genome build.
- Pooled control sample names refer to raw FASTQ sample names, not the pooled output name.
- `homer_genome` is installed in the HOMER environment or points to a FASTA file.
- If HOMER is skipped, check `logs/*.homer.*.log` for the peak count after blacklist filtering.

## Notes for AI Agents

Read `agent/manifest.yml` and `agent/io.contract.yml` before changing workflow logic. Avoid editing raw data, project-specific config files, or sample metadata unless the user explicitly asks. Prefer dry-runs and rule-specific log inspection before modifying rules.
