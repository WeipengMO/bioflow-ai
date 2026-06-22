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

The default workflow writes three browser-track scales:

- `cpm`: CPM-normalized signal. This is useful for routine browser inspection, but it can hide true global R-loop gain or loss because each sample is rescaled to its own mapped-depth total.
- `absolute_spikein`: absolute spike-in calibration with `scale = 10000 / spikein_count(sample)`. This is intended for cross-sample spike-in-normalized display and follows the direct E. coli fragment-count calibration style used in many CUT&Tag workflows.
- `matched_ref_spikein`: matched NoRNaseH-anchor spike-in scaling. Each treatment sample is the anchor for itself and its matched RNase H/control samples, making this scale useful for RNaseH-oriented NoRNaseH versus RNaseH comparisons.

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
scale_methods:
  - cpm
  - absolute_spikein
  - matched_ref_spikein

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

When spike-in is enabled, `absolute_spikein` and `matched_ref_spikein` use the `spikein_count` defined by `spikein.mode`.

`absolute_spikein` and `matched_ref_spikein` intentionally use different constants because they answer different display questions:

```text
absolute_spikein_scale_factor =
  10000 / spikein_fragments(target)

matched_anchor_spikein_scale_factor =
  spikein_fragments(anchor) / spikein_fragments(target)
  * 1000000 / host_unique_fragments(anchor)
```

In this formula, `anchor` means the matched NoRNaseH treatment sample that defines the comparison unit. `target` means the current sample being scaled, which may be the NoRNaseH anchor itself, its matched RNase H sample, or an associated control. If `target == anchor`, the matched-anchor scale reduces to `1000000 / host_unique_fragments(anchor)`. The `1000000` term is a per-million host-fragment display unit for the anchor; it is not the same concept as the `10000` absolute spike-in display constant.

## Sample Metadata Logic

Each FASTQ sample appears once in `config/samples.yml`. The workflow uses the fields below to distinguish R-loop signal samples, RNase H controls, MACS3 background controls, biological replicate groups, and matched pairs.

| Field | Required for | Meaning |
| --- | --- | --- |
| `role` | all samples | One of `treatment`, `rnaseh_control`, or `control`. |
| `group` | `treatment` samples used for replicate-aware outputs | Biological replicate group. Samples with the same `group` are used for consensus/intersect peaks and DESeq2 RNaseH statistics. |
| `control` | optional on `treatment` samples | Matched input/IgG/background sample for MACS3 peak calling. This is not the RNase H comparison sample. |
| `rnaseh_control` | recommended on `treatment` samples | Matched RNase H-treated sample used to define the NoRNaseH/RNaseH pair. |

`role: treatment` means the NoRNaseH/untreated R-loop CUT&Tag sample. `role: rnaseh_control` means the matched RNase H-treated sample. `role: control` means an optional MACS3 background such as input or IgG; it is not used as the RNaseH condition in DESeq2.

The pairing is declared from the `treatment` sample:

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

In this example, `WT_CutTag_RLoop_rep1` is paired with `WT_RNaseH_rep1`. The `WT_RNaseH_rep1` entry itself usually only needs `role: rnaseh_control`; the matched-pair relationship is created by the treatment sample's `rnaseh_control` field.

For replicate-aware DESeq2 RNaseH depletion, each `group` must contain at least two treatment samples with matched `rnaseh_control` samples. The default model is:

```r
design = ~ pair + condition
contrast = NoRNaseH vs RNaseH
```

`pair` blocks each matched treatment/RNaseH pair, and `condition` tests the systematic NoRNaseH versus RNaseH effect. By default, this workflow does not add extra DESeq2 covariates.

## Optional DESeq2 Covariates

Do not add batch-like covariates unless the metadata reflect real experimental design and the covariate is not fully confounded with `condition`. Incorrect covariates can make the model hard to interpret or not full rank.

The current workflow default is no extra covariates. Keep `design_covariates` empty unless the workflow has been extended for covariate-aware DESeq2 designs:

```yaml
rnaseh_sensitive:
  mode: both
  design_covariates: []  # default: no extra DESeq2 covariates
  # design_covariates:
  #   - batch
  #   - prep_batch
```

Optional covariate values, when used in a covariate-aware model, should be short strings without spaces. Use identical spelling for identical batches. If a covariate is listed in `design_covariates`, every sample in the tested NoRNaseH/RNaseH contrast should have a non-empty value for that covariate.

| Covariate | Meaning | Example values | When to consider it |
| --- | --- | --- | --- |
| `batch` | Broad experimental or sequencing batch. | `batch1`, `batch2`, `seq_run_A` | Different experiment dates, operators, flowcells, or sequencing runs may introduce systematic differences. |
| `prep_batch` | Library preparation, CUT&Tag prep, or Tn5/reagent batch. | `prepA`, `prep_2026_06_01`, `tn5_lot_A` | Libraries were prepared in separate prep rounds or with different reagent lots. |
| `cell_fraction` | Biological material or fraction used for the assay. | `nucleus`, `whole_cell`, `chromatin_fraction` | Samples genuinely come from different cell fractions, and each condition has coverage across those fractions. |

Example optional metadata:

```yaml
WT_CutTag_RLoop_rep1:
  role: treatment
  group: WT
  rnaseh_control: WT_RNaseH_rep1
  batch: batch1
  prep_batch: prepA
WT_RNaseH_rep1:
  role: rnaseh_control
  batch: batch1
  prep_batch: prepA
```

Do not encode biological condition as a covariate. For example, do not set `batch: NoRNaseH` for treatment samples and `batch: RNaseH` for RNase H samples. Also avoid adding `batch` if every pair has a unique batch already captured by `pair`, or if all NoRNaseH samples are in one batch and all RNaseH samples are in another batch.

DESeq2 models covariates through the design formula and expects raw integer counts as input; do not manually adjust counts to remove batch effects before DESeq2. See the [Bioconductor DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for design formulas such as `~ batch + condition`, the recommendation to place the variable of interest at the end of the formula, and the warning that DESeq2 input should be un-normalized integer counts. If using DESeq2 in published work, cite Love, Huber, and Anders 2014, Genome Biology, DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).

## Key Outputs

| Output | Meaning |
| --- | --- |
| `results/aligned_data/{sample}.signal.keepdup.filtered.bam` | Filtered duplicate-marked BAM for signal, FRiP, and count quantification |
| `results/aligned_data/{sample}.peaks.dedup.filtered.bam` | Deduplicated/peak-mode BAM when configured for peak calling |
| `results/bigwig/cpm/{sample}.cpm.bw` | CPM-normalized browser signal |
| `results/bigwig/absolute_spikein/{sample}.absolute_spikein.bw` | Absolute spike-in calibrated browser signal using `10000 / spikein_count(sample)` |
| `results/bigwig/matched_ref_spikein/{sample}.matched_ref_spikein.bw` | Matched NoRNaseH-anchor spike-in normalized browser signal |
| `results/qc/normalization/spikein_summary.tsv` | Generic spike-in metrics and scale factors |
| `results/counts/peak_universe.bed` | Shared peak universe |
| `results/counts/peak_universe.saf` | SAF annotation passed to featureCounts |
| `results/counts/peak_counts.featureCounts.txt` | Raw featureCounts output |
| `results/counts/peak_counts.raw.tsv` | Raw shared-universe counts |
| `results/counts/peak_counts.cpm.tsv` | CPM counts |
| `results/counts/peak_counts.spikein_normalized.tsv` | Matched NoRNaseH-anchor spike-in scaled counts, generated only when spike-in is enabled |
| `results/rnaseh/{sample}.rnaseh_sensitive_ratio.bed` | Exploratory RNaseH-sensitive ratio calls |
| `results/rnaseh/{contrast}.rnaseh_depleted.deseq2.tsv` | Replicate-aware RNaseH depletion statistics |
| `results/rnaseh/{contrast}.rnaseh_depleted.bed` | FDR/log2FC-filtered RNaseH-depleted regions |
| `results/qc/replicate_correlation.tsv` | Correlation over shared peak counts |
| `results/qc/replicate_correlation_heatmap.pdf` | Count-correlation heatmap |
| `results/qc/sample_pca.pdf` | PCA over shared peak counts |
| `results/qc/peak_width_distribution.pdf` | Peak width distribution |
| `results/qc/rnaseh_depletion_summary.tsv` | Combined RNaseH ratio/statistical summaries |
| `results/qc/rnaseh_specificity_summary.tsv` | R-loop-specific RNaseH depletion/sensitivity summary |
| `results/qc/rnaseh_signal_scatter.pdf` | No-RNaseH versus RNaseH signal scatter over ratio-tested regions |
| `results/qc/rnaseh_depletion_fraction.pdf` | RNaseH-sensitive/depleted fraction plot |
| `results/qc/spikein_summary.tsv` | QC-level copy/summary of spike-in metrics |
| `results/qc/blacklist_mito_summary.tsv` | Filtering-status summary |
| `results/qc/frip_fragment_level.tsv` | FRiP summary table; paired-end mode uses fragment/read-pair counts |

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
- Host genome Bowtie2 options use `bowtie2_extra`; spike-in Bowtie2 options use the independent `spikein_bowtie2_extra`, which defaults to CUT&Tag-style `--no-overlap --no-dovetail` filtering in paired-end mode.
- Optional TSS/TES/gene-body metaplots require genome annotation and are not enabled by default in this lightweight workflow.
- If spike-in reads are absent or below thresholds, inspect `warnings/cut_tag_rloop_spikein.warning.*`; use `cpm` cautiously as a fallback.
