# CUT&Tag R-loop Profiling Workflow

This BioFlowAI workflow processes CUT&Tag-style R-loop or DNA hybrid-associated libraries from FASTQ files into filtered BAMs, normalized signal tracks, raw sample peaks, shared peak-universe counts, RNaseH-sensitive candidate regions, replicate-aware RNaseH-depleted regions, and R-loop-specific QC summaries.

Use cautious terminology when interpreting the outputs. Raw MACS3 peaks are regions of R-loop-associated CUT&Tag signal; they are not definitive R-loop sites. Stronger R-loop-associated calls should be based on loss of signal after RNase H treatment, ideally quantified over a shared peak universe with replicate-aware statistics. If the assay uses S9.6, antibody specificity can be a major confounder, so RNase H controls and orthogonal validation remain important.

## Biological Output Hierarchy

Recommended interpretation order:

1. `results/peaks/{sample}.cut_tag_rloop_peaks.*Peak`: raw sample peaks for discovery and browser review.
2. `results/intersect_peaks/` and `results/consensus_peaks/`: replicate support summaries.
3. `results/counts/peak_universe.bed`: shared union/consensus peak universe.
4. `results/counts/peak_counts.*.tsv`: common-peak count matrices for quantitative comparison.
5. `results/rnaseh/{sample}.rnaseh_sensitive_ratio.bed`: exploratory RNaseH-sensitive candidates from signal ratios.
6. `results/rnaseh/{contrast}.rnaseh_depleted.deseq2.tsv` and `.bed`: statistically supported RNaseH-depleted regions.
7. Differential R-loop-associated regions should be derived from shared counts and replicate-aware statistics, not sample-specific raw peaks alone.

For global signal changes, spike-in normalization is preferred over CPM because CPM can mask real genome-wide gain or loss of DNA hybrid-associated signal.

## Inputs

Required inputs:

- Paired-end or single-end FASTQs under `raw_data/`, or `config/fastq_manifest.tsv`.
- `config/config.yml`, copied from `config/config.example.yml`.
- `config/samples.yml`, copied from `config/samples.example.yml`.
- Host genome Bowtie2 index prefix and MACS3 genome size.
- Spike-in Bowtie2 index when `spikein.enabled: true`.

Recommended inputs:

- Matched input/background controls for MACS3.
- Matched RNase H-treated controls for each treatment sample.
- Biological replicate groups in `samples.yml`.
- A genome-build matched blacklist BED.

## Configuration Highlights

Main biological config sections:

```yaml
spikein:
  enabled: true
  genome: ecoli
  bowtie2_index: data/spikein/ecoli_mg1655
  mode: proper_pair_fragments
  min_spikein_reads: 1000
  warn_low_fraction: 0.001

count_matrix:
  enabled: true
  peak_universe: consensus
  min_peak_width: 50
  merge_distance: 100
  min_mapq: 30
  featurecounts_extra: ""

rnaseh_sensitive:
  enabled: true
  mode: both
  pseudocount: 0.1
  min_treatment_signal: 0.5
  min_abs_signal_diff: 0.2
  min_fold_change: 2.0
  fdr_threshold: 0.05
  log2fc_threshold: 1.0
```

`spikein.genome` is a label; `spikein.bowtie2_index` controls the actual alignment target. E. coli is only the default example.

Sample metadata keeps the existing role fields:

```yaml
WT_CutTag_RLoop_rep1:
  role: treatment
  group: WT
  control: WT_Input_rep1
  rnaseh_control: WT_RNaseH_rep1
WT_RNaseH_rep1:
  role: rnaseh_control
WT_Input_rep1:
  role: control
```

`group` should identify biological replicate groups used for consensus peaks and replicate-aware RNaseH statistics.

## Key Outputs

| Output | Meaning |
| --- | --- |
| `results/aligned_data/{sample}.signal.keepdup.filtered.bam` | Filtered duplicate-marked BAM for signal, FRiP, and count quantification |
| `results/aligned_data/{sample}.peaks.dedup.filtered.bam` | Deduplicated/peak-mode BAM when configured for peak calling |
| `results/bigwig/CPM/{sample}.CPM.bw` | CPM-normalized browser signal |
| `results/bigwig/spikein/{sample}.spikein.bw` | Spike-in normalized browser signal |
| `results/qc/normalization/spikein_summary.tsv` | Generic spike-in metrics and scale factors |
| `results/counts/peak_universe.bed` | Shared peak universe |
| `results/counts/peak_universe.saf` | SAF annotation passed to featureCounts |
| `results/counts/peak_counts.featureCounts.txt` | Raw featureCounts output |
| `results/counts/peak_counts.raw.tsv` | Raw shared-universe counts |
| `results/counts/peak_counts.cpm.tsv` | CPM counts |
| `results/counts/peak_counts.spikein_normalized.tsv` | Spike-in scaled counts |
| `results/rnaseh/{sample}.rnaseh_sensitive_ratio.bed` | Exploratory RNaseH-sensitive ratio calls |
| `results/rnaseh/{contrast}.rnaseh_depleted.deseq2.tsv` | Replicate-aware RNaseH depletion statistics |
| `results/rnaseh/{contrast}.rnaseh_depleted.bed` | FDR/log2FC-filtered RNaseH-depleted regions |
| `results/qc/replicate_correlation.tsv` | Correlation over shared peak counts |
| `results/qc/replicate_correlation_heatmap.pdf` | Count-correlation heatmap |
| `results/qc/sample_pca.pdf` | PCA over shared peak counts |
| `results/qc/peak_width_distribution.pdf` | Peak width distribution |
| `results/qc/rnaseh_depletion_summary.tsv` | Combined RNaseH ratio/statistical summaries |
| `results/qc/spikein_summary.tsv` | QC-level copy/summary of spike-in metrics |
| `results/qc/blacklist_mito_summary.tsv` | Filtering-status summary |
| `results/qc/frip_fragment_level.tsv` | FRiP summary table |

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/cut_tag_rloop
./scripts/deploy_pipeline.sh /path/to/cut_tag_rloop_project
```

Edit:

```text
/path/to/cut_tag_rloop_project/config/config.yml
/path/to/cut_tag_rloop_project/config/samples.yml
```

Dry-run:

```bash
cd /path/to/cut_tag_rloop_project
./run_snakemake.sh -n -p
```

Run:

```bash
./run_snakemake.sh
```

## Notes And Limitations

- Ratio-based RNaseH-sensitive regions are exploratory and not replicate-aware.
- DESeq2 RNaseH depletion requires at least two treatment/RNaseH pairs in a `group`.
- The count matrix is generated with `featureCounts`; paired-end mode uses `-p --countReadPairs -B -C`.
- Optional TSS/TES/gene-body metaplots require genome annotation and are not enabled by default in this lightweight workflow.
- If spike-in reads are absent or below thresholds, inspect `warnings/cut_tag_rloop_spikein.warning.*`; use CPM cautiously as a fallback.
