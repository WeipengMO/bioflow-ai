# RNA-seq Expression Workflow

## Inputs

Project files:

```text
config/config.yml
config/samples.tsv
raw_data/*.fastq.gz or *.fq.gz
```

`samples.tsv` example:

```tsv
sample	fq1	fq2	group	replicate
control_1	raw_data/control_1_R1.fastq.gz	raw_data/control_1_R2.fastq.gz	control	1
control_2	raw_data/control_2_R1.fastq.gz	raw_data/control_2_R2.fastq.gz	control	2
treated_1	raw_data/treated_1_R1.fastq.gz	raw_data/treated_1_R2.fastq.gz	treated	1
treated_2	raw_data/treated_2_R1.fastq.gz	raw_data/treated_2_R2.fastq.gz	treated	2
```

- `sample`: unique sample name (required)
- `fq1`: R1 FASTQ or single-end FASTQ (required)
- `fq2`: R2 FASTQ for paired-end data; leave empty for single-end data (optional)
- `group`: sample group for downstream analysis (optional)
- `replicate`: replicate ID for distinguishing biological or technical replicates (optional, values can be 1, 2, 3, etc.)
- Extra columns are allowed for downstream use

If `group` or `replicate` columns are not provided, the workflow will still run correctly. These columns are only used for downstream analysis or metadata purposes.

### Reference Resource Preparation

All genome-related files must use the same assembly.

### STAR Index

`genome.star_index` should be a STAR-built index directory. For example:

```yaml
genome:
  gtf: path/to/hg38.gtf
  star_index: path/to/hg38_star_index
```

Build the STAR index:

```bash
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir path/to/hg38_star_index \
  --genomeFastaFiles path/to/hg38.fa \
  --sjdbGTFfile path/to/hg38.gtf
```

### GTF Annotation

`genome.gtf` should be a standard GTF file.

## Configuration

Main config options (config/config.yml):

```yaml
genome:
  gtf: /path/to/annotation.gtf
  star_index: /path/to/star_index

strandness: unstranded  # unstranded / forward / reverse
```

- `strandness` controls featureCounts `-s` parameter:
  - `unstranded` → `-s 0`
  - `forward` → `-s 1`
  - `reverse` → `-s 2`

## Outputs

```text
results/matrix/gene_counts.tsv      # Raw count matrix, with gene_symbol column
results/matrix/gene_fpkm.tsv        # FPKM matrix, with gene_symbol column
results/matrix/gene_tpm.tsv         # TPM matrix, with gene_symbol column
results/bam/{sample}.sorted.bam     # Aligned BAM
results/bam/{sample}.sorted.bam.bai
results/bigwig/{sample}.bw          # bigWig tracks
results/report/multiqc_report.html  # QC summary
```

## Workflow Logic

1. `fastp`: Raw FASTQ QC and filtering
2. `align_star`: STAR alignment, outputs BAM
3. `count_featurecounts`: featureCounts quantification, generates raw counts
4. `matrix`: Merge counts, compute FPKM/TPM, output matrices with gene_symbol column
5. `bigwig`: bamCoverage to generate bigWig
6. `report`: MultiQC summary report

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/rnaseq
./scripts/deploy_pipeline.sh /path/to/rnaseq_project
cd /path/to/rnaseq_project
```

Place FASTQ files under `raw_data/`, edit `config/config.yml` and `config/samples.tsv`, then run:

```bash
./run_snakemake.sh -n
./run_snakemake.sh
```

The wrapper uses Snakemake conda environments by default and installs them under:

```text
/data/user/mowp/snakemake_conda_envs
```

Set it persistently in `config/run_snakemake.env` if needed:

```bash
SNAKEMAKE_CONDA_PREFIX=/path/to/snakemake_conda_envs
```

`config/run_snakemake.env.example` is only a template and is not loaded by `run_snakemake.sh`.

Or override this per run:

```bash
SNAKEMAKE_CONDA_PREFIX=/path/to/snakemake_conda_envs ./run_snakemake.sh
```

## Troubleshooting

- Ensure sample names match between FASTQ files and samples.tsv
- STAR index must be pre-built and match the GTF
- gene_counts.tsv, gene_fpkm.tsv, and gene_tpm.tsv all include a gene_symbol column
- Always start with a dry-run:

```bash
./run_snakemake.sh -n -p
```

## Notes for AI Agents

- Read this README and the files under `agent/` before modifying config or starting a run.
- Prefer deploying the workflow into a separate project directory instead of editing files under `workflows/rnaseq/` directly.
- Use `scripts/deploy_pipeline.sh` or the documented deployment steps when creating a new analysis directory.
- Preserve existing `raw_data/` files and treat user-provided sample tables and config files as user-owned inputs unless the user explicitly requests edits.
- Start with `./run_snakemake.sh -n -p` before attempting a full run.
- When reporting back, summarize created files, unresolved required inputs, and the next recommended command.


## External Resource Links

- STAR: https://github.com/alexdobin/STAR
- Subread/featureCounts: http://subread.sourceforge.net/
- deepTools: https://deeptools.readthedocs.io/
- MultiQC: https://multiqc.info/

For downstream analysis such as differential expression or enrichment, use tools like DESeq2, edgeR, or clusterProfiler.