# ATAC-seq Workflow

This BioFlowAI workflow processes ATAC-seq FASTQ files into filtered BAM files, normalized bigWig tracks, MACS3 narrow peaks, replicate consensus peaks, and core ATAC-seq QC metrics.

The default logic follows a mainstream bulk ATAC-seq path:

1. FASTQ quality control and adapter trimming with `fastp`.
2. Alignment to a Bowtie2 genome index.
3. Coordinate sorting and duplicate removal.
4. Post-alignment filtering by MAPQ, SAM flags, and mitochondrial chromosomes.
5. Duplicate removal.
6. Optional blacklist-region removal.
7. CPM-normalized bigWig generation.
8. MACS3 narrow peak calling.
9. Union-style replicate consensus peaks for groups with at least two biological replicates.
10. Insert-size, FRiP, optional TSS enrichment, and MultiQC reporting.

## Inputs

Required inputs:

- Paired-end or single-end FASTQ files under `raw_data/`, or a manifest at `config/fastq_manifest.tsv`.
- `config/config.yml`, copied from `config/config.example.yml`.
- `config/samples.yml`, copied from `config/samples.example.yml`.
- A Bowtie2 genome index prefix.
- MACS3 genome size, such as `hs`, `mm`, or an effective genome size number.

Recommended inputs:

- ENCODE blacklist BED for the same genome build.
- TSS BED for TSS enrichment profiling.

FASTQ auto-discovery supports common paired-end names such as:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
sample_1.fq.gz
sample_2.fq.gz
```

For irregular names or unmerged lanes, generate and edit a manifest:

```bash
python scripts/make_fastq_manifest.py --raw-data-dir raw_data --output config/fastq_manifest.tsv
```

Then set:

```yaml
fastq_manifest: config/fastq_manifest.tsv
```

## Configuration

Main config keys:

| Key | Meaning |
| --- | --- |
| `mode` | `pe` or `se`; paired-end is recommended for ATAC-seq |
| `raw_data_dir` | FASTQ directory used when no manifest is provided |
| `sample_config` | Sample metadata YAML |
| `genome` | Bowtie2 index prefix |
| `gsize` | MACS3 genome size |
| `blacklist` | Optional blacklist BED removed after duplicate marking |
| `mitochondrial_chromosomes` | Chromosome names removed before duplicate marking |
| `tss_bed` | TSS BED used by deepTools when TSS enrichment is enabled |
| `enable_tss_enrichment` | Produce per-sample TSS enrichment profiles |
| `enable_multiqc` | Produce a MultiQC report |
| `fastp_extra` | Extra fastp arguments used by both PE and SE |
| `fastp_pe_extra` | PE-only fastp arguments |
| `fastp_se_extra` | SE-only fastp arguments |

Sample metadata uses one mapping entry per sample:

```yaml
WT_ATAC_rep1:
  group: WT_ATAC
WT_ATAC_rep2:
  group: WT_ATAC
```

`group` is optional. Groups with at least two samples produce a union-style consensus peak set at `results/replicate_intersect/{group}_intersect.bed`.

## Outputs

Default output root: `results/`.

Important outputs:

| Output | Description |
| --- | --- |
| `aligned_data/*.sorted.rmdup.filtered.bam` | Deduplicated, filtered BAM files |
| `bigwig/*.CPM.bw` | CPM-normalized signal tracks |
| `macs3_results/narrow/*_peaks.narrowPeak` | MACS3 narrow peaks |
| `replicate_intersect/*_intersect.bed` | Union-style consensus peaks retained by replicate support |
| `qc/fastp/` | FASTQ QC reports |
| `qc/mark_duplicates/` | Duplicate metrics |
| `qc/insert_size/` | ATAC fragment-size metrics and histograms for paired-end mode |
| `qc/frip/*.frip.txt` | FRiP summary tables |
| `qc/tss_enrichment/` | Optional TSS profile outputs |
| `reports/multiqc_report.html` | MultiQC summary |
| `logs/` | Rule logs |

## Workflow Logic

The analysis BAM is built as:

```text
FASTQ
  -> fastp
  -> bowtie2 | samtools sort
  -> mapped/proper-pair/MAPQ/mitochondrial chromosome filtering
  -> Picard MarkDuplicates
  -> optional blacklist removal
```

Default alignment filters before duplicate removal:

- `alignment_filter_min_mapq: 30`
- paired-end: `samtools view -f 2 -F 3852`
- single-end: `samtools view -F 3844`

Default MACS3 behavior:

- paired-end: `-f BAMPE --keep-dup all -q 0.01`
- single-end: `--nomodel --shift -100 --extsize 200 --keep-dup all -q 0.01`

Duplicates are already removed before peak calling, so `--keep-dup all` prevents MACS3 from applying an additional duplicate model.

Replicate consensus logic:

1. Clean each replicate peak file to BED3 and merge within each replicate.
2. Concatenate all replicate peak files.
3. Merge them into a candidate peak universe.
4. Count how many replicate peak sets support each candidate peak.
5. Keep candidates with support `>= consensus_min_support`.

This produces a union-style consensus peak set that is suitable as a shared feature universe for downstream differential accessibility analysis.

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/atacseq
./scripts/deploy_pipeline.sh /path/to/atacseq_project
```

Then edit:

```text
/path/to/atacseq_project/config/config.yml
/path/to/atacseq_project/config/samples.yml
```

Run a dry-run first:

```bash
cd /path/to/atacseq_project
./run_snakemake.sh -n -p
```

Run the workflow:

```bash
./run_snakemake.sh
```

The wrapper reads `config/run_snakemake.env`, which must define:

```bash
SNAKEMAKE_CORES=16
SNAKEMAKE_CONDA_PREFIX=/path/to/snakemake_conda_envs
```

## Troubleshooting

- Missing FASTQs: check sample names in `config/samples.yml` against FASTQ names or `config/fastq_manifest.tsv`.
- Empty filtered BAMs: check `alignment_filter_min_mapq`, `alignment_filter_view_extra`, mitochondrial chromosome names, and blacklist genome build.
- No TSS output: set `enable_tss_enrichment: true` and provide `tss_bed`.
- MACS3 empty peaks: confirm `gsize`, read depth, filtering stringency, and whether the library has expected ATAC insert-size periodicity.

## Notes for AI Agents

- Read this README and `agent/` before deploying or changing logic.
- Prefer `scripts/deploy_pipeline.sh` for project setup.
- Do not edit source templates under `workflows/atacseq/` when the task is only to run a project-specific analysis.
- Always dry-run before a full run.
- Keep genome index, blacklist, TSS BED, and `gsize` on the same genome build.
