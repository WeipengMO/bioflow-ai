# CUT&Tag R-loop Profiling Workflow

This BioFlowAI workflow processes CUT&Tag-style R-loop profiling libraries from FASTQ files into filtered signal BAMs, normalized bigWig tracks, MACS3 CUT&Tag R-loop peaks, RNaseH-sensitive regions, replicate peak sets, normalization metrics, and QC reports.

The default logic is intentionally scoped to CUT&Tag R-loop analysis:

1. FASTQ quality control and adapter trimming with `fastp`.
2. Human genome alignment with `bowtie2`.
3. MAPQ, SAM flag, mitochondrial chromosome, and optional blacklist filtering.
4. Duplicate marking for QC while keeping duplicate-marked fragments for signal and peak calling.
5. CPM and E. coli spike-in normalized bigWig generation under method-specific `bigwig/` subdirectories.
6. E. coli spike-in alignment and group-local spike-in scale factor calculation.
7. MACS3 peak calling with `--keep-dup all`.
8. RNaseH-sensitive region identification from treatment peaks by treatment-vs-RNaseH signal depletion.
9. Strict replicate `intersect_peaks` and support-based `consensus_peaks`.
10. FRiP, insert-size metrics, normalization metrics, spike-in warnings, and MultiQC reporting.

## Inputs

Required inputs:

- Paired-end or single-end FASTQ files under `raw_data/`, or an explicit manifest at `config/fastq_manifest.tsv`.
- `config/config.yml`, copied from `config/config.example.yml`.
- `config/samples.yml`, copied from `config/samples.example.yml`.
- A human Bowtie2 genome index prefix.
- MACS3 genome size, such as `hs`, `mm`, or a numeric effective genome size.
- E. coli spike-in reads in the same FASTQs and an E. coli Bowtie2 index prefix for default spike-in scaling.

Recommended inputs:

- Matched input/background controls for MACS3 when available.
- Matched RNase H-treated libraries for RNaseH-sensitive region identification.
- ENCODE blacklist BED for the same human genome build.

FASTQ auto-discovery supports common names such as:

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
| `assay` | Workflow label; default `CUT&Tag-R-loop` |
| `mode` | `pe` or `se` |
| `raw_data_dir` | FASTQ directory used when no manifest is provided |
| `sample_config` | Sample metadata YAML |
| `genome` | Human Bowtie2 index prefix |
| `gsize` | MACS3 genome size |
| `blacklist` | Optional blacklist BED removed from BAMs |
| `peak_type` | `broad` or `narrow`; broad is the default for R-loop-enriched domains |
| `scale_methods` | List containing `CPM`, `spikein`, or both; default is both |
| `spikein_genome` | E. coli Bowtie2 index prefix; required by the default workflow |
| `spikein_min_mapped_reads` | Warning threshold for low E. coli mapped reads |
| `spikein_min_fraction` | Warning threshold for low E. coli read fraction |
| `peak_duplicate_mode` | `auto`, `keep_marked`, or `remove`; `auto` keeps duplicate-marked BAMs for CUT&Tag peak calling |
| `rnaseh_signal_min_fold_change` | Minimum treatment/RNaseH signal ratio for RNaseH-sensitive regions; default `2.0` |
| `rnaseh_signal_min_treatment_signal` | Minimum treatment signal for RNaseH-sensitive regions; default `0.0` |
| `consensus_min_support` | Minimum replicate support for consensus peak sets |

Sample metadata uses one mapping entry per FASTQ sample:

```yaml
WT_CutTag_RLoop_rep1:
  role: treatment
  group: WT_CutTag_RLoop_rep1
  control: WT_Input_rep1
  rnaseh_control: WT_RNaseH_rep1
WT_Input_rep1:
  role: control
WT_RNaseH_rep1:
  role: rnaseh_control
```

- `role: treatment` marks CUT&Tag R-loop samples that should be peak-called.
- `role: control` marks input/background controls used by MACS3.
- `role: rnaseh_control` marks matched RNase H-treated libraries.
- `group` names the spike-in normalization group for one treatment/RNaseH pair.
- `control` points from a treatment sample to a MACS3 control sample.
- `rnaseh_control` points from a treatment sample to its matched RNase H-treated sample.

## Scale Methods

`scale_methods` is a list. The default workflow runs both supported methods:

- `CPM`: depth-normalized signal tracks.
- `spikein`: group-local E. coli spike-in scaled signal tracks.

Example:

```yaml
scale_methods:
  - CPM
  - spikein
spikein_genome: data/spikein/ecoli_mg1655
```

For spike-in scaling, each treatment sample's `group` defines one normalization group and the treatment sample is the reference for its matched RNase H control. The final bigWig scale factor is:

```text
spikein_reads_reference / spikein_reads_sample * 1e6 / human_unique_fragments_reference
```

With the default `spikein` method, the workflow writes:

```text
results/qc/normalization/normalization_metrics.tsv
warnings/cut_tag_rloop_spikein.warning.tsv
warnings/cut_tag_rloop_spikein.warning.txt
```

`normalization_metrics.tsv` records:

```text
sample
human_total_reads
human_mapped_reads
human_unique_fragments
ecoli_total_reads
ecoli_mapped_reads
ecoli_unique_fragments
ecoli_fraction
spikein_group
spikein_reference_sample
spikein_reference_ecoli_mapped_reads
spikein_reference_human_unique_fragments
spikein_raw_scale_factor
spikein_unit_scale_factor
spikein_scale_factor
scale_method
warning_level
warning_message
```

If any sample has `ecoli_mapped_reads=0`, the spike-in branch fails after writing warning files outside `results/`. Low nonzero spike-in support writes warnings and continues.

## Outputs

Default output root: `results/`.

Important outputs:

| Output | Description |
| --- | --- |
| `aligned_data/{sample}.signal.keepdup.filtered.bam` | Filtered duplicate-marked signal BAM |
| `bigwig/CPM/{sample}.CPM.bw` | CPM-normalized bigWig |
| `bigwig/spikein/{sample}.spikein.bw` | Group-local spike-in normalized bigWig |
| `peaks/{sample}.cut_tag_rloop_peaks.*Peak` | MACS3 CUT&Tag R-loop peaks |
| `rnaseh_sensitive/{method}/{sample}.rnaseh_sensitive_signal.tsv` | Per-peak treatment/RNaseH signal table |
| `rnaseh_sensitive/{method}/{sample}.rnaseh_sensitive_regions.bed` | RNaseH-sensitive regions passing signal thresholds |
| `intersect_peaks/{group}.intersect_peaks.bed` | Strict all-replicate supported peaks |
| `consensus_peaks/{group}.consensus_peaks.bed` | Support-threshold consensus peaks |
| `rnaseh_sensitive_consensus/{method}/{group}.rnaseh_sensitive_consensus.bed` | Replicate-level RNaseH-sensitive consensus regions |
| `qc/normalization/normalization_metrics.tsv` | Human/E. coli normalization metrics and spike-in scale factors |
| `qc/frip/{sample}.frip.txt` | FRiP over MACS3 CUT&Tag R-loop peaks |
| `qc/rnaseh_sensitive/{method}/{sample}.summary.tsv` | RNaseH-sensitive region summary |
| `reports/multiqc_report.html` | MultiQC summary |
| `warnings/cut_tag_rloop_spikein.warning.*` | Project-level spike-in warning files outside `results/` |

The bigWig files are signal tracks produced by different normalization methods. Downstream analysis should center on peak files, RNaseH-sensitive regions, consensus peak sets, and counts over consensus peaks rather than treating bigWig tracks as the primary quantitative result.

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/cut_tag_rloop
./scripts/deploy_pipeline.sh /path/to/cut_tag_rloop_project
```

Then edit:

```text
/path/to/cut_tag_rloop_project/config/config.yml
/path/to/cut_tag_rloop_project/config/samples.yml
```

Run a dry-run first:

```bash
cd /path/to/cut_tag_rloop_project
./run_snakemake.sh -n -p
```

Run the workflow:

```bash
./run_snakemake.sh
```

## Troubleshooting

- Missing FASTQs: check sample names in `config/samples.yml` against FASTQ names or `config/fastq_manifest.tsv`.
- Missing spike-in index: build one with `bash scripts/build_spikein_index.sh`, then set `spikein_genome`.
- Empty filtered BAMs: check genome build, `alignment_filter_min_mapq`, SAM flag filters, mitochondrial chromosome names, and blacklist file.
- Weak signal: inspect input/RNase H controls, duplicate rate, library complexity, and peak mode.
- RNaseH-sensitive output from CPM should be treated as exploratory when global depletion is the biological question.
- Spike-in warnings are written outside `results/` in `warnings/`; check them before interpreting spike-in scaled results.
