# R-loop Profiling Workflow

This BioFlowAI workflow processes R-loop profiling libraries such as DRIP-seq, DRIPc-seq, R-ChIP, MapR, and related S9.6 or RNase H based assays from FASTQ files into filtered BAM files, normalized bigWig tracks, MACS3 peak/domain calls, replicate consensus peaks, optional RNase H sensitivity-filtered peaks, and QC summaries.

The default logic follows a mainstream bulk R-loop analysis path:

1. FASTQ quality control and adapter trimming with `fastp`.
2. Alignment to a Bowtie2 genome index.
3. MAPQ, SAM flag, and optional mitochondrial chromosome filtering.
4. Duplicate marking/removal.
5. Optional blacklist-region removal.
6. CPM-normalized bigWig generation.
7. MACS3 broad peak/domain calling by default, with optional narrow peak mode.
8. Optional RNase H control subtraction to retain R-loop-sensitive peaks.
9. Union-style replicate consensus peaks by biological group.
10. FRiP, optional fragment-size metrics, and MultiQC reporting.

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
| `enable_rnaseh_subtraction` | Remove peaks overlapping sample-specific RNase H controls when metadata provides `rnaseh_control` |
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

## Outputs

Default output root: `results/`.

Important outputs:

| Output | Description |
| --- | --- |
| `aligned_data/*.sorted.rmdup.filtered.bam` | Deduplicated, filtered BAM files for all samples |
| `bigwig/*.CPM.bw` | CPM-normalized signal tracks |
| `macs3_results/{broad,narrow}/*_peaks.*Peak` | MACS3 R-loop domains or peaks for treatment samples |
| `rnaseh_filtered/*_rnaseh_sensitive.*Peak` | Optional RNase H-sensitive treatment peaks |
| `replicate_consensus/*_consensus.bed` | Union-style treatment-group consensus peaks |
| `qc/fastp/` | FASTQ QC reports |
| `qc/mark_duplicates/` | Duplicate metrics |
| `qc/frip/*.frip.txt` | FRiP summary tables for treatment samples |
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

- broad mode: `macs3 callpeak --broad --broad-cutoff 0.1 -q 0.05`
- narrow mode: `macs3 callpeak -q 0.01`
- paired-end uses `-f BAMPE`; single-end uses `-f BAM`

For DRIP-seq and related broad-domain assays, keep `peak_type: broad` unless the experimental design specifically expects sharp sites.

Replicate consensus logic:

1. Use RNase H-sensitive peaks when `enable_rnaseh_subtraction: true` and the treatment has `rnaseh_control` metadata; otherwise use the primary MACS3 peak.
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
- Empty RNase H-sensitive peaks: confirm that treatment and RNase H control libraries are comparable and use the same genome build.
- Always start with a dry-run.

## Notes for AI Agents

- Read this README and `agent/` before deploying or changing logic.
- Prefer `scripts/deploy_pipeline.sh` for project setup.
- Do not edit source templates under `workflows/rloop/` when the task is only to run a project-specific analysis.
- Treat raw data, sample metadata, and project config as user-owned inputs.
- Keep Bowtie2 index, blacklist, and MACS3 `gsize` on the same genome build.

## External Resource Links

- MACS3: https://macs3-project.github.io/MACS/
- Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/
- deepTools: https://deeptools.readthedocs.io/
- MultiQC: https://multiqc.info/
