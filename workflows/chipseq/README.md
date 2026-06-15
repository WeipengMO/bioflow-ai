# ChIP-seq / CUT&Tag Workflow

## Inputs

Project files:

```text
config/config.yml
config/samples.yml
raw_data/*.fastq.gz or raw_data/*.fq.gz
```

Optional files:

```text
config/fastq_manifest.tsv
config/chipqc_sample_sheet.tsv
```

Reference files configured in `config/config.yml`:

| Key | Required when | Meaning |
|---|---:|---|
| `genome` | always | Bowtie2 index prefix passed to `bowtie2 -x` |
| `gsize` | peak calling | MACS3 genome size, such as `hs`, `mm`, or a numeric effective genome size |
| `gtf` | `enable_homer: true` | Gene annotation GTF passed to `annotatePeaks.pl -gtf` |
| `picard_path` | always | `picard` command or `/path/to/picard.jar` for MarkDuplicates |
| `regions_bed` | `enable_deeptools_profile: true` | BED regions used by `computeMatrix` |
| `homer_genome` | `enable_homer: true` | HOMER genome name such as `hg38`, or a FASTA path accepted by HOMER |
| `blacklist` | optional blacklist filtering | BED file reused for final BAM filtering and HOMER peak filtering |

### Reference Files

All genome-related files must use the same genome assembly. Do not mix `hg19`, `hg38`, `mm10`, and other assemblies in the same run.

### Bowtie2 Index

`genome` must be a Bowtie2 index prefix, not only a FASTA path. For example, if these files exist:

```text
/path/to/hg38/bowtie2_index/hg38.1.bt2
/path/to/hg38/bowtie2_index/hg38.2.bt2
...
```

set:

```yaml
genome: /path/to/hg38/bowtie2_index/hg38
```

You can build the index from a genome FASTA:

```bash
bowtie2-build /path/to/hg38.fa /path/to/hg38/bowtie2_index/hg38
```

### regions_bed for deepTools

`regions_bed` is the region set used by:

```bash
computeMatrix ... -R ${regions_bed} -S sample.bw
```

Common choices are gene bodies, promoters, enhancers, or a custom BED of loci relevant to the experiment. For a UCSC RefSeq gene-body BED on hg38:

```bash
wget -O /path/to/hg38.refGene.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat /path/to/hg38.refGene.txt.gz \
  | awk 'BEGIN{OFS="\t"} {print $3,$5,$6,$13,0,$4}' \
  | sort -k1,1 -k2,2n \
  > /path/to/hg38.refGene.bed
```

Then set:

```yaml
regions_bed: /path/to/hg38.refGene.bed
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
homer_genome: /path/to/hg38.fa
```

### GTF for HOMER Annotation

`gtf` is used by the HOMER peak annotation step:

```bash
annotatePeaks.pl peaks.bed ${homer_genome} -gtf ${gtf}
```

Set it to a gene annotation GTF from the same assembly as `genome`, `regions_bed`, `blacklist`, and `homer_genome`:

```yaml
gtf: /path/to/annotations.gtf
```

### Blacklist

`blacklist` is optional. If set, both BAM and peak filtering use:

```bash
bedtools intersect -v -a peaks.bed -b ${blacklist}
```

For human hg38, use the ENCODE blacklist from the Boyle-Lab Blacklist repository:

```bash
wget -O /path/to/hg38-blacklist.v2.bed.gz \
  https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip -c /path/to/hg38-blacklist.v2.bed.gz > /path/to/hg38-blacklist.v2.bed
```

Then set:

```yaml
blacklist: /path/to/hg38-blacklist.v2.bed
```

For other assemblies, check the `lists/` directory in the Boyle-Lab Blacklist repository and use the matching assembly file.

### Sample Metadata

`config/samples.yml` defines sample metadata. `config/config.yml` points to it with `sample_config`.

Rules for project setup:

- Each FASTQ sample is listed once.
- `role` is required and must be `treatment` or `control`.
- `role: treatment` marks an experimental ChIP/CUT&Tag sample. These samples are used for MACS3 peak calling, intersect peaks, consensus peaks, ChIPQC rows, and HOMER targets.
- `role: control` marks an input/IgG/control sample. These samples are preprocessed into BAM/bigWig files, but they are not peak-called as treatments.
- `group` names biological replicates for consensus peaks.
- `control` on a treatment points to either one control sample or one control pool.
- `pool` on a control sample is optional. Only control samples with the same `pool` value are merged into a pooled control.

**Example pooled control setup:**

```yaml
WT_H3K27ac_rep1:
  role: treatment
  group: WT_H3K27ac
  control: WT_input_pool
WT_H3K27ac_rep2:
  role: treatment
  group: WT_H3K27ac
  control: WT_input_pool
WT_Input_rep1:
  role: control
  pool: WT_input_pool
WT_Input_rep2:
  role: control
  pool: WT_input_pool
```

This means:

```text
WT_input_pool = merge(WT_Input_rep1, WT_Input_rep2)
WT_H3K27ac_rep1 uses WT_input_pool
WT_H3K27ac_rep2 uses WT_input_pool
```

Here, `WT_H3K27ac_rep1` and `WT_H3K27ac_rep2` are treatment samples and will be peak-called. `WT_Input_rep1` and `WT_Input_rep2` are control samples and are only used through `WT_input_pool`.

**Example matched control setup:**

```yaml
WT_H3K27ac_rep1:
  role: treatment
  group: WT_H3K27ac
  control: WT_Input_rep1
WT_H3K27ac_rep2:
  role: treatment
  group: WT_H3K27ac
  control: WT_Input_rep2
WT_Input_rep1:
  role: control
WT_Input_rep2:
  role: control
```

No `pool` field is needed for matched controls; the treatment points directly to one control sample.
In this example, only the two `role: treatment` samples are peak-called; the two `role: control` samples are used as MACS3 controls.

**Example without controls:**

```yaml
WT_H3K27ac_rep1:
  role: treatment
  group: WT_H3K27ac
WT_H3K27ac_rep2:
  role: treatment
  group: WT_H3K27ac
KO_H3K27ac_rep1:
  role: treatment
  group: KO_H3K27ac
KO_H3K27ac_rep2:
  role: treatment
  group: KO_H3K27ac
```

For a no-control run, do not add `role: control` samples and do not set `control` on treatment samples. In `config/config.yml`, use no-control peak modes:

```yaml
call_peak_modes:
  - without_control
peak_set_mode: without_control
chipqc_peak_mode: without_control
homer_peak_mode: without_control
```

No-control MACS3 outputs are written under `<outdir>/macs3_results/narrow_no_control/` or `<outdir>/macs3_results/broad_no_control/`.

Groups with fewer than two valid treatment samples are ignored.

### FASTQ Input

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

## Configuration

| Key | Meaning |
|---|---|
| `mode` | `pe` or `se` |
| `outdir` | Root directory for generated outputs; default is `results` |
| `sample_config` | YAML file defining sample roles, replicate groups, and controls |
| `raw_data_dir` | Directory containing FASTQ files |
| `fastq_manifest` | Optional TSV with `sample`, `read1`, and optional `read2` columns |
| `genome` | Bowtie2 index prefix |
| `regions_bed` | BED regions for deepTools profile plots |
| `gtf` | Gene annotation GTF used by HOMER peak annotation |
| `gsize` | MACS3 genome size |
| `picard_path` | Picard command or JAR path |
| `call_peak_modes` | `with_control`, `without_control`, or both |
| `call_peak_types` | `narrow`, `broad`, or both |
| `enable_alignment_filter` | Apply MAPQ/flag filtering before `MarkDuplicates`; default `true` |
| `alignment_filter_min_mapq` | Minimum MAPQ passed to `samtools view -q`; default `30` |
| `alignment_filter_view_extra` | Extra arguments appended to pre-MarkDuplicates `samtools view`; if empty, the workflow picks mode-specific defaults |
| `enable_blacklist_filter` | Remove blacklist overlaps after `MarkDuplicates`; default `true` |
| `macs3_broad_cutoff` | Value passed to `--broad-cutoff` for broad peak calling |
| `macs3_extra` | Extra arguments appended directly to `macs3 callpeak`, such as `--qvalue 0.01` |
| `peak_set_mode` | Peak mode used by `intersect_peaks` and `consensus_peaks` |
| `peak_set_type` | Peak type used by `intersect_peaks` and `consensus_peaks` |
| `consensus_min_support_count` | Minimum number of treatment samples required for a consensus peak; default `1` |
| `consensus_min_support_fraction` | Minimum treatment-sample fraction required for a consensus peak; default `0.0` |
| `enable_deeptools_profile` | Generate deepTools matrix/profile plots |
| `enable_chipqc` | Generate ChIPQC report |
| `chipqc_sample_sheet` | Optional user-provided ChIPQC sample sheet |
| `chipqc_peak_mode` | Peak mode used by ChIPQC |
| `chipqc_peak_type` | Peak type used by ChIPQC |
| `enable_homer` | Run HOMER motif enrichment |
| `homer_input_source` | `intersect_peaks` or `sample_peaks` |
| `homer_genome` | HOMER genome name or FASTA path |
| `blacklist` | Optional BED blacklist reused by final BAM filtering and HOMER |
| `homer_min_peaks` | Minimum peaks required to run HOMER; default is 50 |
| `threads` | Per-rule thread settings |

Advanced optional keys supported by the workflow but omitted from `config.example.yml`: `chipqc_annotation`, `chipqc_report_dir`, `chipqc_report_facet`, `homer_report_dir`, `homer_input_dir`, `homer_min_peaks`, and `workdir`.

Default pre-MarkDuplicates alignment filtering is equivalent to:

```bash
samtools view -q 30 -f 2 -F 1804 -b input.sorted.bam > input.sorted.filtered.bam
```

In single-end mode, the default becomes `samtools view -q 30 -F 1796 -b`, which removes unmapped, secondary, QC-fail, and duplicate reads without applying a proper-pair filter.

When `enable_blacklist_filter: true`, duplicate-removed BAMs are filtered again with:

```bash
bedtools intersect -v -abam input.sorted.rmdup.bam -b ${blacklist} > input.sorted.rmdup.filtered.bam
```

## Outputs

```text
<outdir>/clean_data/                         # temporary cleaned FASTQs
<outdir>/aligned_data/*.sorted.rmdup.bam     # duplicate-removed BAMs
<outdir>/aligned_data/*.sorted.rmdup.bam.bai
<outdir>/aligned_data/*.sorted.rmdup.filtered.bam   # final analysis BAMs when enable_blacklist_filter: true
<outdir>/aligned_data/*.sorted.rmdup.filtered.bam.bai
<outdir>/qc/fastp/
<outdir>/qc/mark_duplicates/
<outdir>/bigwig/*.CPM.bw
<outdir>/deeptools_profile/*.scale.png
<outdir>/macs3_results/
<outdir>/intersect_peaks/
<outdir>/consensus_peaks/
<outdir>/reports/chipqc/
<outdir>/reports/homer/
<outdir>/logs/
```

## Workflow Logic

The workflow DAG is built from `config/config.yml`, `config/samples.yml`, and discovered FASTQs or `config/fastq_manifest.tsv`.

Execution order:

1. `fastp`
   - Inputs: raw FASTQs from `raw_data_dir` or `fastq_manifest`.
  - Outputs: temporary cleaned FASTQs under `<outdir>/clean_data/`.
  - Reports: `<outdir>/qc/fastp/{sample}.html` and `.json`.

2. `align_reads`
   - Runs `bowtie2` against `genome`.
   - Adds read group tags (`ID`, `SM`, and `PL`) during alignment.
   - Sorts alignments with `samtools sort`.
  - Outputs temporary sorted BAMs under `<outdir>/aligned_data/`.

3. `filter_aligned`
  - Runs by default before duplicate removal.
  - Applies `samtools view` filtering with `alignment_filter_min_mapq` and `alignment_filter_view_extra`.
  - Can be disabled with `enable_alignment_filter: false`.
  - Outputs temporary alignment-filtered BAMs:

```text
<outdir>/aligned_data/{sample}.sorted.filtered.bam
<outdir>/aligned_data/{sample}.sorted.filtered.bam.bai
```

4. `mark_duplicates`
   - Runs Picard `MarkDuplicates`.
  - Uses the alignment-filtered BAM when `enable_alignment_filter: true`, otherwise uses the sorted BAM.
  - Outputs duplicate-removed BAMs and indexes:

```text
<outdir>/aligned_data/{sample}.sorted.rmdup.bam
<outdir>/aligned_data/{sample}.sorted.rmdup.bam.bai
<outdir>/qc/mark_duplicates/{sample}.metrics.txt
```

5. `filter_blacklist`
  - Runs by default after `mark_duplicates`.
  - If `blacklist` is set, removes blacklist overlaps with `bedtools intersect -v -abam`.
  - Can be disabled with `enable_blacklist_filter: false`.
  - Outputs:

```text
<outdir>/aligned_data/{sample}.sorted.rmdup.filtered.bam
<outdir>/aligned_data/{sample}.sorted.rmdup.filtered.bam.bai
```

6. `bam_coverage`
   - Runs deepTools `bamCoverage`.
  - Uses the blacklist-filtered BAM when `enable_blacklist_filter: true`, otherwise uses the duplicate-removed BAM.
   - Outputs normalized bigWig tracks:

```text
<outdir>/bigwig/{sample}.sorted.rmdup.CPM.bw
```

7. `compute_matrix_profile`
   - Runs only when `enable_deeptools_profile: true`.
   - Uses `regions_bed`.
   - Outputs:

```text
<outdir>/deeptools_profile/{sample}.scale.png
```

8. `merge_pooled_control`
   - Runs only when control samples declare a shared `pool`.
   - Runs once for each pool name.
  - Merges the downstream analysis BAMs for control samples in that pool into:

```text
<outdir>/aligned_data/{pool_name}.sorted.rmdup.bam
<outdir>/aligned_data/{pool_name}.sorted.rmdup.bam.bai
```

When `enable_blacklist_filter: true`, pooled controls are written as:

```text
<outdir>/aligned_data/{pool_name}.sorted.rmdup.filtered.bam
<outdir>/aligned_data/{pool_name}.sorted.rmdup.filtered.bam.bai
```

9. `macs3_callpeak`
   - Runs one job for each configured treatment, peak mode, and peak type.
   - With-control peaks require `matched` or `pooled` controls.
   - No-control peaks use `without_control` and run MACS3 without a `-c` control BAM.
  - Uses the blacklist-filtered BAM when `enable_blacklist_filter: true`.
   - `macs3_extra` is appended directly to `macs3 callpeak` for advanced options.
   - Output directories:

```text
<outdir>/macs3_results/narrow/
<outdir>/macs3_results/broad/
<outdir>/macs3_results/narrow_no_control/
<outdir>/macs3_results/broad_no_control/
```

10. `intersect_peaks`
   - Uses treatment samples that share the same `group`.
   - Intersects replicate peak files with `bedtools multiinter`.
   - Runs only for groups with at least two treatment samples.
   - Consensus intervals must be supported by all replicate peak files in that group.
   - This is a strict all-replicate intersection for conservative QC, not the default downstream peak universe.
   - Outputs:

```text
<outdir>/intersect_peaks/{group}_intersect_peaks.bed
```

11. `consensus_peaks_group` and `consensus_peaks_all_treatments`
   - Use the MACS3 peak files selected by `peak_set_mode` and `peak_set_type`.
   - Merge overlapping or directly adjacent peak intervals into a consensus peak universe.
   - Group-level outputs use treatment samples that share the same `group`.
   - All-treatment outputs use every `role: treatment` sample and exclude controls.
   - `consensus_min_support_count` and `consensus_min_support_fraction` can filter weakly supported union peaks.
   - BED outputs contain merged intervals with `consensus_peak_N` IDs and support counts.
   - TSV outputs record `peak_id`, coordinates, `support_count`, `support_fraction`, and `support_samples`.
   - Support matrix outputs record each treatment sample as present/absent for each consensus peak.
   - SAF outputs can be used by counting tools such as featureCounts.
   - Outputs:

```text
<outdir>/consensus_peaks/groups/{group}_consensus_peaks.bed
<outdir>/consensus_peaks/groups/{group}_consensus_peaks.tsv
<outdir>/consensus_peaks/groups/{group}_consensus_support_matrix.tsv
<outdir>/consensus_peaks/groups/{group}_consensus_peaks.saf
<outdir>/consensus_peaks/all_treatments_consensus_peaks.bed
<outdir>/consensus_peaks/all_treatments_consensus_peaks.tsv
<outdir>/consensus_peaks/all_treatments_consensus_support_matrix.tsv
<outdir>/consensus_peaks/all_treatments_consensus_peaks.saf
```

12. `count_consensus_peaks`
   - Counts final analysis BAM reads over the all-treatment consensus peak universe with `bedtools multicov`.
   - Includes treatment and control samples as count columns, so the matrix can feed downstream differential binding workflows.
   - Outputs:

```text
<outdir>/consensus_peaks/counts/all_treatments_counts.tsv
```

13. `generate_chipqc_sample_sheet`
   - Runs only when `enable_chipqc: true` and `chipqc_sample_sheet` is empty.
   - Derives the sheet from `samples`, control mapping, treatment groups, and selected peak files.
   - Generates:

```text
<outdir>/reports/chipqc/chipqc_sample_sheet.generated.tsv
```

14. `chipqc_report`
    - Runs only when `enable_chipqc: true`.
    - Uses treatment BAMs, control BAMs if configured, selected MACS3 peak files, and the ChIPQC sample sheet.
  - Uses the blacklist-filtered BAM when `enable_blacklist_filter: true`.
    - Outputs the ChIPQC report directory and completion marker:

```text
<outdir>/reports/chipqc/
<outdir>/reports/chipqc/.chipqc_complete
```

15. `homer_prepare_motif_peaks`
    - Runs only when `enable_homer: true`.
    - Uses either intersect peaks or per-sample MACS3 peaks depending on `homer_input_source`.
    - Uses MACS summits for narrowPeak motif analysis when `homer_use_summit: true`; consensus BED inputs fall back to interval centers.
    - If `blacklist` is set, removes blacklisted intervals before motif enrichment.
    - Outputs prepared temporary BED files under `<outdir>/homer_inputs/`.

16. `homer_prepare_annotation_peaks`
    - Runs only when `enable_homer: true`.
    - Uses full peak intervals instead of summit-centered 1 bp intervals.
    - If `blacklist` is set, removes blacklisted intervals before annotation.
    - Outputs prepared temporary BED files under `<outdir>/homer_inputs/`.

17. `homer_find_motifs`
    - Runs `findMotifsGenome.pl`.
    - Skips HOMER if fewer than `homer_min_peaks` peaks remain after filtering.
    - Starts from `homer_prepare_motif_peaks`.
    - Shares the same target report directory with `homer_annotate_peaks`, so the workflow only removes stale motif-specific outputs at rule start.
    - Outputs:

```text
<outdir>/reports/homer/{homer_input_source}_{homer_peak_mode}_{homer_peak_type}/{target}/.motifs_complete
```

18. `homer_annotate_peaks`
    - Runs `annotatePeaks.pl` on the full-interval prepared BED from `homer_prepare_annotation_peaks`.
    - Independent from `homer_find_motifs`; it does not wait for motif discovery outputs.
    - Always passes `-gtf {gtf}` from `config/config.yml`.
    - Writes `<outdir>/reports/homer/.../{target}/annotatePeaks.txt`.

19. `homer_plot_annotation_distribution`
    - Reads `annotatePeaks.txt` and collapses HOMER annotations into major feature classes.
    - Writes a pie chart for promoter, exon, intron, intergenic, TTS, and other categories.
    - Writes `<outdir>/reports/homer/.../{target}/peak_feature_distribution.pie.png`.
    - Also writes `<outdir>/reports/homer/.../{target}/peak_feature_distribution.tsv` with counts and fractions.

17. `homer_plot_annotation_distribution_summary`
    - Collects all per-target `peak_feature_distribution.tsv` files into one table.
    - Writes `<outdir>/reports/homer/.../peak_feature_distribution.summary.tsv`.
    - Writes `<outdir>/reports/homer/.../peak_feature_distribution.stacked_bar.png` with figure height scaled to the number of targets.

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/chipseq
./scripts/deploy_pipeline.sh /path/to/chipseq_project
cd /path/to/chipseq_project
```

Put FASTQ files under `raw_data/`, edit `config/config.yml` and `config/samples.yml`, then run:

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

Always start with a dry-run:

```bash
./run_snakemake.sh -n -p
```

Common checks:

- Sample names in `config/samples.yml` match FASTQ-derived names or the `sample` column in `config/fastq_manifest.tsv`.
- `genome` points to a Bowtie2 index prefix, not only a `.fa` file.
- `regions_bed`, `blacklist`, `homer_genome`, and `genome` use the same assembly.
- For pooled controls, treatment `control` values point to a pool name, and control samples declare that same `pool`.
- For no-control runs, treatment samples must omit `control`, and `call_peak_modes` plus downstream peak modes should use `without_control`.
- `call_peak_modes` includes `with_control` only when a matched or pooled control is configured.
- ChIPQC uses the selected `chipqc_peak_mode` and `chipqc_peak_type`; those files must also be produced by `call_peak_modes` and `call_peak_types`.
- If HOMER is skipped, check `<outdir>/logs/*.homer.prepare.*.log` and `<outdir>/logs/*.homer.*.log` for the number of peaks after blacklist filtering.

## Notes for AI Agents

- Read this README and the files under `agent/` before editing config or running commands.
- Prefer deploying the workflow into a separate project directory instead of modifying files under `workflows/chipseq/` directly.
- Use `scripts/deploy_pipeline.sh` or the documented deployment steps when preparing a new analysis directory.
- Treat `raw_data/`, user config files, and sample metadata as user-owned inputs unless the user explicitly asks for changes.
- Start with `./run_snakemake.sh -n -p` before attempting a real execution.
- When reporting back, summarize created files, unresolved required inputs, and the next recommended command.

## External Resource Links

- ENCODE blacklist files: <https://github.com/Boyle-Lab/Blacklist>
- hg38 blacklist v2 BED gzip: <https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz>
- UCSC downloads FAQ: <https://genome.ucsc.edu/FAQ/FAQdownloads.html>
- UCSC hg38 RefSeq table: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz>
