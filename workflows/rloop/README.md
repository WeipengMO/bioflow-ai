# R-loop Profiling Workflow

This BioFlowAI workflow processes R-loop profiling libraries such as DRIP-seq, DRIPc-seq, R-ChIP, MapR, and related S9.6 or RNase H based assays from FASTQ files into filtered BAM files, CPM browser bigWig tracks, optional common-scale signal tracks, MACS3 peak/domain calls, replicate consensus peaks, RNase H no-overlap filtered peaks, quantitative RNaseH signal summaries, and QC reports.

The default logic follows a mainstream bulk R-loop analysis path:

1. FASTQ quality control and adapter trimming with `fastp`.
2. Alignment to a Bowtie2 genome index.
3. MAPQ, SAM flag, and optional mitochondrial chromosome filtering.
4. Duplicate marking/removal.
5. Optional blacklist-region removal.
6. CPM-normalized bigWig generation for local browser shape comparisons.
7. Optional common-scale bigWig generation for global signal comparisons such as RNaseH depletion.
8. MACS3 broad peak/domain calling by default, with optional narrow peak mode.
9. Optional RNase H no-overlap filtering with `bedtools intersect -v`.
10. Quantitative treatment-vs-RNaseH signal summaries over no-overlap peaks.
11. Union-style replicate consensus peaks by biological group.
12. FRiP, optional fragment-size metrics, bigWig header QC, and MultiQC reporting.

## Inputs

Required inputs:

- Paired-end or single-end FASTQ files under `raw_data/`, or an explicit manifest at `config/fastq_manifest.tsv`.
- `config/config.yml`, copied from `config/config.example.yml`.
- `config/samples.yml`, copied from `config/samples.example.yml`.
- A Bowtie2 genome index prefix.
- MACS3 genome size, such as `hs`, `mm`, or a numeric effective genome size.

Recommended inputs:

- Matched input, IgG, RNase H-treated, or catalytically inactive RNase H control libraries when available.
- ENCODE blacklist BED for the same genome build.
- GTF annotation for downstream interpretation outside this workflow.

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
| `assay` | Free-text assay label, for example `DRIP-seq`, `DRIPc-seq`, `R-ChIP`, or `MapR` |
| `mode` | `pe` or `se` |
| `raw_data_dir` | FASTQ directory used when no manifest is provided |
| `sample_config` | Sample metadata YAML |
| `genome` | Bowtie2 index prefix |
| `gsize` | MACS3 genome size |
| `peak_type` | `broad` or `narrow`; broad is recommended for DRIP-seq domains |
| `blacklist` | Optional blacklist BED removed after duplicate marking and from peaks |
| `enable_rnaseh_subtraction` | Remove treatment peaks overlapping sample-specific RNase H control peaks with `bedtools intersect -v`; this is no-overlap filtering, not statistical depletion testing |
| `enable_common_scale_bigwig` | Also produce common-scale bigWigs for between-sample signal comparisons |
| `common_scale_normalization` | `None`/`raw` for `bamCoverage --normalizeUsing None --scaleFactor`, or `RPGC` for depth-normalized coverage with an explicit scale factor |
| `signal_scale_factors_tsv` | TSV with `sample` and `scale_factor`; required when `enable_common_scale_bigwig: true` |
| `signal_scale_factor_method` | `effective_fragments`, `raw`, `tsv`, or `spikein` |
| `spikein_genome` | Bowtie2 index prefix for the spike-in reference genome |
| `spikein_method` | `ratio` or `rpm`; how spike-in counts become scale factors |
| `spikein_reference_sample` | Reference sample for ratio mode (gets scale_factor=1.0) |
| `peak_duplicate_mode` | `auto` (default), `keep_marked`, or `remove`; auto uses markdup BAM with `--keep-dup all` |
| `dup_threshold_for_dedup` | Warn when Picard PERCENT_DUPLICATION >= this (default 0.90); auto mode only |
| `rnaseh_signal_min_fold_change` | Minimum treatment/RNaseH signal ratio for the depleted BED output; default `2.0` |
| `rnaseh_signal_min_treatment_signal` | Minimum treatment signal for the depleted BED output; default `0.0` |
| `enable_insert_size_qc` | Produce Picard insert-size metrics in paired-end mode |
| `enable_multiqc` | Produce a MultiQC report |

Sample metadata uses one mapping entry per FASTQ sample:

```yaml
WT_DRIP_rep1:
  role: treatment
  group: WT_DRIP
  control: WT_Input_rep1
  rnaseh_control: WT_RNaseH_rep1
WT_Input_rep1:
  role: control
WT_RNaseH_rep1:
  role: rnaseh_control
```

- `role: treatment` marks samples that should be peak-called.
- `role: control` marks input, IgG, or other background controls used by MACS3.
- `role: rnaseh_control` marks RNase H-treated control libraries used for optional specificity subtraction.
- `group` names biological replicate groups for consensus peaks.
- `control` points from a treatment sample to a MACS3 control sample.
- `rnaseh_control` points from a treatment sample to the matching RNase H control sample.

### Signal Track Modes

The workflow intentionally separates browser visualization from RNaseH-depletion interpretation:

| Track mode | Output | Intended use |
| --- | --- | --- |
| CPM browser track | `results/bigwig/*.CPM.bw` | Local shape comparison in a genome browser |
| Common-scale track | `results/bigwig_common_scale/*.scaled.bw` | Between-sample signal comparison where global depletion must be preserved |

Per-sample CPM rescales each library independently. That is appropriate for seeing local enrichment shape, but it can hide a global RNaseH signal drop by forcing RNaseH and untreated samples onto comparable total signal. For RNaseH-depletion interpretation, provide spike-in, input, or author-defined scale factors in `signal_scale_factors_tsv` and enable common-scale bigWigs.

Example scale-factor file:

```text
sample	scale_factor	normalization_group	note
RLoop_rep1	1.0	RLoop	reference scale
RNaseH_rep1	0.18	RLoop	spike-in or processed-track-derived scale
```

### Spike-in Normalization

For experiments where global signal differences must be preserved (e.g., RNaseH depletion analysis), spike-in normalization aligns reads to an E. coli K-12 MG1655 reference genome and computes per-sample scale factors from spike-in alignment counts.

Setup:

```bash
bash scripts/build_spikein_index.sh  # downloads E. coli MG1655 and builds Bowtie2 index under data/spikein/
```

Config:

```yaml
signal_scale_factor_method: spikein
spikein_genome: data/spikein/ecoli_mg1655
spikein_method: ratio               # or rpm
spikein_reference_sample: WT_DRIP_rep1  # required for ratio mode
enable_common_scale_bigwig: true
```

Two modes:

- **ratio**: `scale_factor = reference_spikein_count / sample_spikein_count`. The reference sample gets `scale_factor = 1.0`. Samples with fewer spike-in reads get a larger scale factor, amplifying their signal to match the reference.
- **rpm**: `scale_factor = 1,000,000 / sample_spikein_count`. Normalizes to reads per million spike-in reads.

### Duplicate Handling for CUT&Tag/CUT&RUN

CUT&Tag and CUT&Run PCR duplicates often enrich in real signal regions. Removing them can bias toward background noise. The default `peak_duplicate_mode: auto` uses the duplicate-marked (but not removed) BAM for peak calling with MACS3 `--keep-dup all`, while still running Picard MarkDuplicates for QC metrics.

Only consider `peak_duplicate_mode: remove` when the Picard duplication rate is extremely high (80-90%+), indicating near sequencing saturation.

## Outputs

Default output root: `results/`.

Important outputs:

| Output | Description |
| --- | --- |
| `aligned_data/*.sorted.rmdup.filtered.bam` | Deduplicated, filtered BAM files for all samples |
| `bigwig/*.CPM.bw` | CPM-normalized browser tracks for local shape inspection |
| `bigwig_common_scale/*.scaled.bw` | Optional common-scale signal tracks for between-sample comparisons |
| `macs3_results/{broad,narrow}/*_peaks.*Peak` | MACS3 R-loop domains or peaks for treatment samples |
| `rnaseh_no_overlap/*_rnaseh_no_overlap.*Peak` | Treatment peaks with RNase H control peak overlaps removed by `bedtools intersect -v` |
| `rnaseh_filtered/*_rnaseh_sensitive.*Peak` | Deprecated compatibility alias for `rnaseh_no_overlap`; this file does not by itself prove RNaseH sensitivity |
| `rnaseh_signal/*_rnaseh_signal.tsv` | Per-peak treatment and RNaseH signal, ratio, difference, and classification |
| `rnaseh_signal/*_rnaseh_depleted.bed` | Peaks passing configurable treatment > RNaseH and fold-change thresholds |
| `replicate_consensus/*_consensus.bed` | Union-style treatment-group consensus peaks |
| `qc/fastp/` | FASTQ QC reports |
| `qc/mark_duplicates/` | Duplicate metrics |
| `qc/frip/*.frip.txt` | FRiP summary tables for treatment samples |
| `qc/bigwig_scale/*.bigwig_header.tsv` | bigWig header summaries including `sumData`, covered bases, max value, and scale factor |
| `qc/rnaseh_sensitivity/*.summary.tsv` | Raw peak count, no-overlap retained count, overlap-removed fraction, and signal-depleted count |
| `reports/multiqc_report.html` | MultiQC summary |
| `logs/` | Rule logs |

## Workflow Logic

The analysis BAM is built as:

```text
FASTQ
  -> fastp
  -> bowtie2 | samtools sort
  -> MAPQ / SAM flag / optional mitochondrial filtering
  -> Picard MarkDuplicates
  -> optional blacklist removal
```

Default filtering:

- `alignment_filter_min_mapq: 30`
- paired-end: `samtools view -F 3852`
- single-end: `samtools view -F 3844`

Default R-loop peak calling:

- broad mode: `macs3 callpeak --broad --broad-cutoff 0.1 --keep-dup all -q 0.05`
- narrow mode: `macs3 callpeak --keep-dup all -q 0.01`
- paired-end uses `-f BAMPE`; single-end uses `-f BAM`
- `--keep-dup all` is always passed explicitly regardless of `peak_duplicate_mode`
- Default `peak_duplicate_mode: auto` uses markdup BAM (duplicates retained) for peak calling

For DRIP-seq and related broad-domain assays, keep `peak_type: broad` unless the experimental design specifically expects sharp sites.

RNaseH-aware peak logic:

1. When `enable_rnaseh_subtraction: true` and a treatment has `rnaseh_control` metadata, the final peak set uses `rnaseh_no_overlap`, meaning treatment peaks that do not overlap merged RNase H control peaks. This is overlap filtering only.
2. The deprecated `rnaseh_filtered/*_rnaseh_sensitive.*Peak` path is written as a compatibility alias to the same no-overlap file.
3. Quantitative RNaseH depletion is assessed separately in `rnaseh_signal/*_rnaseh_signal.tsv` and `rnaseh_signal/*_rnaseh_depleted.bed`.

Replicate consensus logic:

1. Use RNase H no-overlap peaks when `enable_rnaseh_subtraction: true` and the treatment has `rnaseh_control` metadata; otherwise use the primary MACS3 peak.
2. Clean each replicate peak to BED3 and merge within each replicate.
3. Build a merged candidate peak universe.
4. Count replicate support for each candidate.
5. Keep candidates with support `>= consensus_min_support`.

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/rloop
./scripts/deploy_pipeline.sh /path/to/rloop_project
```

Then edit:

```text
/path/to/rloop_project/config/config.yml
/path/to/rloop_project/config/samples.yml
```

Run a dry-run first:

```bash
cd /path/to/rloop_project
./run_snakemake.sh -n -p
```

Run the workflow:

```bash
./run_snakemake.sh
```

## Troubleshooting

- Missing FASTQs: check sample names in `config/samples.yml` against FASTQ names or `config/fastq_manifest.tsv`.
- Empty filtered BAMs: check genome build, `alignment_filter_min_mapq`, SAM flag filters, mitochondrial chromosome names, and blacklist file.
- Weak R-loop signal: inspect input/RNase H controls, duplicate rate, library complexity, and whether the assay should be modeled as broad or narrow.
- Empty RNase H no-overlap peaks: confirm that treatment and RNase H control libraries are comparable, peak files use the same genome build, and the RNase H control is not behaving like an untreated positive sample.
- RNaseH appears similar to treatment in CPM bigWigs: CPM rescales each sample independently. Enable common-scale bigWigs with spike-in or supplied scale factors before interpreting global depletion.
- Weak spike-in signal: verify spike-in genome matches the organism used in library prep, check Bowtie2 alignment rate in logs, ensure spike-in DNA was added before library prep.
- High duplication rate warning: check Picard markdup metrics; if >90%, consider `peak_duplicate_mode: remove` as a sensitivity analysis.
- Always start with a dry-run.

## Notes for AI Agents

- Read this README and `agent/` before deploying or changing logic.
- Prefer `scripts/deploy_pipeline.sh` for project setup.
- Do not edit source templates under `workflows/rloop/` when the task is only to run a project-specific analysis.
- Treat raw data, sample metadata, and project config as user-owned inputs.
- Keep Bowtie2 index, blacklist, and MACS3 `gsize` on the same genome build.

## External Resource Links

- MACS3 callpeak: https://macs3-project.github.io/MACS/docs/callpeak.html
- Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/
- deepTools bamCoverage: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
- deepTools bamCompare: https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html
- MultiQC: https://multiqc.info/
