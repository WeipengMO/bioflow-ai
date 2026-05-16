# ChIP-seq / CUT&Tag Workflow

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/chipseq
./scripts/deploy_pipeline.sh /path/to/chipseq_project
cd /path/to/chipseq_project
```

Prepare project config and metadata:

```bash
cp config/config.example.yml config/config.yml
cp config/samples.example.yml config/samples.yml
cp config/replicates.example.yml config/replicates.yml
mkdir -p raw_data resources
```

Put FASTQ files under `raw_data/`, edit `config/config.yml`, `config/samples.yml`, and `config/replicates.yml`, then run:

```bash
./run_snakemake.sh -n
./run_snakemake.sh
```

The wrapper uses Snakemake conda environments by default and installs them under:

```text
/data/user/mowp/snakemake_conda_envs
```

Override this per run if needed:

```bash
SNAKEMAKE_CONDA_PREFIX=/path/to/snakemake_conda_envs ./run_snakemake.sh
```

## Required Inputs

Project files:

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

Reference resources configured in `config/config.yml`:

| Key | Required when | Meaning |
|---|---:|---|
| `genome` | always | Bowtie2 index prefix passed to `bowtie2 -x` |
| `gsize` | peak calling | MACS3 genome size, such as `hs`, `mm`, or a numeric effective genome size |
| `picard_path` | always | `picard` command or `/path/to/picard.jar` for MarkDuplicates |
| `regions_bed` | `enable_deeptools_profile: true` | BED regions used by `computeMatrix` |
| `homer_genome` | `enable_homer: true` | HOMER genome name such as `hg38`, or a FASTA path accepted by HOMER |
| `homer_blacklist` | optional HOMER filtering | BED file of regions excluded before motif enrichment |
| `chipqc_annotation` | optional ChIPQC annotation | Annotation value passed to `ChIPQC()`; leave empty to omit |

## Reference Resource Preparation

All genome-related files must use the same genome assembly. Do not mix `hg19`, `hg38`, `mm10`, and other assemblies in the same run.

### Bowtie2 Index

`genome` must be a Bowtie2 index prefix, not only a FASTA path. For example, if these files exist:

```text
resources/hg38/bowtie2_index/hg38.1.bt2
resources/hg38/bowtie2_index/hg38.2.bt2
...
```

set:

```yaml
genome: resources/hg38/bowtie2_index/hg38
```

You can build the index from a genome FASTA:

```bash
bowtie2-build resources/hg38.fa resources/hg38/bowtie2_index/hg38
```

### regions_bed for deepTools

`regions_bed` is the region set used by:

```bash
computeMatrix ... -R ${regions_bed} -S sample.bw
```

Common choices are gene bodies, promoters, enhancers, or a custom BED of loci relevant to the experiment. For a UCSC RefSeq gene-body BED on hg38:

```bash
mkdir -p resources
wget -O resources/hg38.refGene.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat resources/hg38.refGene.txt.gz \
  | awk 'BEGIN{OFS="\t"} {print $3,$5,$6,$13,0,$4}' \
  | sort -k1,1 -k2,2n \
  > resources/hg38.refGene.bed
```

Then set:

```yaml
regions_bed: resources/hg38.refGene.bed
```

For other assemblies, replace `hg38` in the UCSC URL with the target assembly if the table is available.

### MACS3 Genome Size

`gsize` is passed to `macs3 callpeak -g`. Typical values:

```yaml
gsize: hs   # human
gsize: mm   # mouse
```

You can also provide a numeric effective genome size if your project standard requires it.

### HOMER Genome

`homer_genome` is passed to:

```bash
findMotifsGenome.pl peaks.bed ${homer_genome} outdir
```

Use a HOMER-installed genome name, for example:

```yaml
homer_genome: hg38
```

or a FASTA path if that is how your HOMER environment is configured:

```yaml
homer_genome: resources/hg38.fa
```

### HOMER Blacklist

`homer_blacklist` is optional. If set, peaks are filtered with:

```bash
bedtools intersect -v -a peaks.bed -b ${homer_blacklist}
```

For human hg38, use the ENCODE blacklist from the Boyle-Lab Blacklist repository:

```bash
mkdir -p resources
wget -O resources/hg38-blacklist.v2.bed.gz \
  https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip -c resources/hg38-blacklist.v2.bed.gz > resources/hg38-blacklist.v2.bed
```

Then set:

```yaml
homer_blacklist: resources/hg38-blacklist.v2.bed
```

For other assemblies, check the `lists/` directory in the Boyle-Lab Blacklist repository and use the matching assembly file.

## Sample Metadata

`config/samples.yml` defines which samples are treatments and how controls are used.

Supported control strategies:

| Strategy | Meaning |
|---|---|
| `none` | No control BAM is used; only no-control peaks are valid |
| `matched` | Each treatment maps to one control in `controls.matched` |
| `pooled` | One or more control groups are merged, and each treatment is assigned to exactly one pooled control |

Example pooled control setup:

```yaml
treatments:
  - WT_H3K27ac_rep1
  - WT_H3K27ac_rep2
  - KO_H3K27ac_rep1
  - KO_H3K27ac_rep2

controls:
  strategy: pooled

  pooled:
    WT_input_pool:
      treatments:
        - WT_H3K27ac_rep1
        - WT_H3K27ac_rep2
      controls:
        - WT_Input_rep1
        - WT_Input_rep2

    KO_input_pool:
      treatments:
        - KO_H3K27ac_rep1
        - KO_H3K27ac_rep2
      controls:
        - KO_Input_rep1
        - KO_Input_rep2
```

This means:

```text
WT_input_pool = merge(WT_Input_rep1, WT_Input_rep2)
WT_H3K27ac_rep1 uses WT_input_pool
WT_H3K27ac_rep2 uses WT_input_pool

KO_input_pool = merge(KO_Input_rep1, KO_Input_rep2)
KO_H3K27ac_rep1 uses KO_input_pool
KO_H3K27ac_rep2 uses KO_input_pool
```

Example matched control setup:

```yaml
treatments:
  - WT_H3K27ac_rep1
  - WT_H3K27ac_rep2
  - KO_H3K27ac_rep1
  - KO_H3K27ac_rep2

controls:
  strategy: matched
  matched:
    WT_H3K27ac_rep1: WT_Input_rep1
    WT_H3K27ac_rep2: WT_Input_rep2
    KO_H3K27ac_rep1: KO_Input_rep1
    KO_H3K27ac_rep2: KO_Input_rep2
```

`config/replicates.yml` defines biological replicate groups used for consensus peaks:

```yaml
WT_H3K27ac:
  - WT_H3K27ac_rep1
  - WT_H3K27ac_rep2
```

Groups with fewer than two valid treatment samples are ignored.

## FASTQ Input

For paired-end data, automatic FASTQ discovery supports names such as:

```text
sample_R1.fastq.gz / sample_R2.fastq.gz
sample_1.fq.gz     / sample_2.fq.gz
sample.R1.fq.gz    / sample.R2.fq.gz
sample-R1.fq.gz    / sample-R2.fq.gz
sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
```

For single-end data, set:

```yaml
mode: se
```

For non-standard names or unmerged multi-lane FASTQs, provide a manifest:

```tsv
sample	read1	read2
WT_rep1	raw_data/a_L001_R1.fastq.gz	raw_data/a_L001_R2.fastq.gz
```

Generate a draft manifest:

```bash
python scripts/make_fastq_manifest.py \
  --raw-data-dir raw_data \
  --output config/fastq_manifest.tsv
```

Then set:

```yaml
fastq_manifest: config/fastq_manifest.tsv
```

## Running Logic

The workflow DAG is built from `config/config.yml`, `config/samples.yml`, discovered FASTQs or `config/fastq_manifest.tsv`, and optional replicate groups.

Execution order:

1. `fastp`
   - Inputs: raw FASTQs from `raw_data_dir` or `fastq_manifest`.
   - Outputs: temporary cleaned FASTQs under `clean_data/`.
   - Reports: `qc/fastp/{sample}.html` and `.json`.

2. `align_reads`
   - Runs `bowtie2` against `genome`.
   - Sorts alignments with `samtools sort`.
   - Outputs temporary sorted BAMs under `aligned_data/`.

3. `mark_duplicates`
   - Runs Picard `MarkDuplicates`.
   - Outputs final duplicate-removed BAMs and indexes:

```text
aligned_data/{sample}.sorted.rmdup.bam
aligned_data/{sample}.sorted.rmdup.bam.bai
qc/mark_duplicates/{sample}.metrics.txt
```

4. `bam_coverage`
   - Runs deepTools `bamCoverage`.
   - Outputs normalized bigWig tracks:

```text
tracks/{sample}.sorted.rmdup.CPM.bw
```

5. `compute_matrix_profile`
   - Runs only when `enable_deeptools_profile: true`.
   - Uses `regions_bed`.
   - Outputs:

```text
deeptools_profile/{sample}.scale.png
```

6. `merge_pooled_control`
   - Runs only when `controls.strategy: pooled`.
   - Runs once for each entry under `controls.pooled`.
   - Merges `controls.pooled.<pool_name>.controls` into:

```text
aligned_data/{pool_name}.sorted.rmdup.bam
aligned_data/{pool_name}.sorted.rmdup.bam.bai
```

7. `macs3_callpeak`
   - Runs one job for each configured treatment, peak mode, and peak type.
   - With-control peaks require `matched` or `pooled` controls.
   - Output directories:

```text
macs3_results/narrow/
macs3_results/broad/
macs3_results/narrow_no_control/
macs3_results/broad_no_control/
```

8. `replicate_intersect`
   - Uses `config/replicates.yml`.
   - Intersects replicate peak files with `bedtools multiinter`.
   - `replicate_min_support` controls the minimum number of replicates required; if omitted, all replicates in the group must support the interval.
   - Outputs:

```text
replicate_intersect/{group}_intersect.bed
```

9. `generate_chipqc_sample_sheet`
   - Runs only when `enable_chipqc: true` and `chipqc_sample_sheet` is empty.
   - Generates:

```text
reports/chipqc/chipqc_sample_sheet.generated.tsv
```

10. `chipqc_report`
    - Runs only when `enable_chipqc: true`.
    - Uses treatment BAMs, control BAMs if configured, selected MACS3 peak files, and the ChIPQC sample sheet.
    - Outputs the ChIPQC report directory and completion marker:

```text
reports/chipqc/
reports/chipqc/.chipqc_complete
```

11. `homer_prepare_peaks`
    - Runs only when `enable_homer: true`.
    - Uses either replicate consensus peaks or per-sample MACS3 peaks depending on `homer_input_source`.
    - If `homer_blacklist` is set, removes blacklisted intervals before motif enrichment.
    - Outputs prepared temporary BED files under `homer_inputs/`.

12. `homer_find_motifs`
    - Runs `findMotifsGenome.pl`.
    - Skips HOMER if fewer than `homer_min_peaks` peaks remain after filtering.
    - Outputs:

```text
reports/homer/{homer_input_source}_{homer_peak_mode}_{homer_peak_type}/{target}/
```

## Main Outputs

```text
clean_data/                         # temporary cleaned FASTQs
aligned_data/*.sorted.rmdup.bam     # final BAMs
aligned_data/*.sorted.rmdup.bam.bai
qc/fastp/
qc/mark_duplicates/
tracks/*.CPM.bw
deeptools_profile/*.scale.png
macs3_results/
replicate_intersect/
reports/chipqc/
reports/homer/
logs/
```

## Configuration Keys

| Key | Meaning |
|---|---|
| `mode` | `pe` or `se` |
| `raw_data_dir` | Directory containing FASTQ files |
| `fastq_manifest` | Optional TSV with `sample`, `read1`, and optional `read2` columns |
| `samples` | YAML file defining treatment and control samples |
| `replicates` | YAML file defining replicate groups |
| `genome` | Bowtie2 index prefix |
| `regions_bed` | BED regions for deepTools profile plots |
| `gsize` | MACS3 genome size |
| `picard_path` | Picard command or JAR path |
| `call_peak_modes` | `with_control`, `without_control`, or both |
| `call_peak_types` | `narrow`, `broad`, or both |
| `replicate_peak_mode` | Peak mode used for replicate consensus |
| `replicate_peak_type` | Peak type used for replicate consensus |
| `replicate_min_support` | Optional minimum replicate support for consensus intervals |
| `enable_deeptools_profile` | Generate deepTools matrix/profile plots |
| `enable_chipqc` | Generate ChIPQC report |
| `chipqc_sample_sheet` | Optional user-provided ChIPQC sample sheet |
| `chipqc_peak_mode` | Peak mode used by ChIPQC |
| `chipqc_peak_type` | Peak type used by ChIPQC |
| `chipqc_annotation` | Optional annotation argument for ChIPQC |
| `enable_homer` | Run HOMER motif enrichment |
| `homer_input_source` | `replicate_intersect` or `sample_peaks` |
| `homer_genome` | HOMER genome name or FASTA path |
| `homer_blacklist` | Optional BED blacklist before HOMER |
| `homer_min_peaks` | Minimum peaks required to run HOMER; default is 50 |
| `threads` | Per-rule thread settings |

## Troubleshooting

Always start with a dry-run:

```bash
./run_snakemake.sh -n -p
```

Common checks:

- Sample names in `config/samples.yml` match FASTQ-derived names or the `sample` column in `config/fastq_manifest.tsv`.
- `genome` points to a Bowtie2 index prefix, not only a `.fa` file.
- `regions_bed`, `homer_blacklist`, `homer_genome`, and `genome` use the same assembly.
- In pooled mode, every treatment must appear in exactly one `controls.pooled.<pool_name>.treatments` list.
- `controls.pooled.<pool_name>.controls` contains raw control sample names, not the pooled output name.
- `call_peak_modes` includes `with_control` only when a matched or pooled control is configured.
- ChIPQC uses the selected `chipqc_peak_mode` and `chipqc_peak_type`; those files must also be produced by `call_peak_modes` and `call_peak_types`.
- If HOMER is skipped, check `logs/*.homer.prepare.*.log` and `logs/*.homer.*.log` for the number of peaks after blacklist filtering.

## External Resource Links

- ENCODE blacklist files: <https://github.com/Boyle-Lab/Blacklist>
- hg38 blacklist v2 BED gzip: <https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz>
- UCSC downloads FAQ: <https://genome.ucsc.edu/FAQ/FAQdownloads.html>
- UCSC hg38 RefSeq table: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz>
