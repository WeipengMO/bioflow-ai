# BioFlowAI Hi-C / HiChIP workflow

This workflow provides an agent-ready Snakemake template for Hi-C and HiChIP analysis. It is designed to be deployed into a project directory with its reusable workflow files linked from `workflows/hic_hichip/`.

## Supported assays

- **Hi-C**: optional FASTQ trimming/QC, HiC-Pro processing, contact matrices, `.cool/.mcool`, QC, loop calling through Mustache or cooltools, and differential loop analysis.
- **HiChIP**: optional FASTQ trimming/QC, HiC-Pro processing, FitHiChIP-oriented loop calling, loop count matrix generation, differential loop analysis, and optional promoter/enhancer annotation.

## Workflow logic

Main steps:

1. Optional `fastp` paired-end read trimming and QC. This is disabled by default; enable it only when FASTQ QC shows adapter or low-quality tail issues.
2. HiC-Pro processing through Singularity.
3. Link HiC-Pro `allValidPairs` and matrix outputs into standardized BioFlowAI result paths.
4. Convert HiC-Pro sparse matrices to `.cool` and `.mcool`.
5. Collect HiC-Pro QC and generate summary plots/MultiQC report.
6. Call loops per sample.
7. Build a consensus loop universe.
8. Quantify loop contacts across samples.
9. Run differential loop analysis using DESeq2 when biological replicates are available; otherwise write descriptive fold-change results.
10. Optionally annotate promoter/enhancer loops.

## Required inputs

Edit `config/config.yml` and provide:

- Paired-end FASTQ files.
- `config/samples.tsv` with at least: `sample`, `fq1`, `fq2`, `assay`, `condition`, `replicate`, `group`.
- `config/comparisons.tsv` for differential loop analysis.
- Reference genome FASTA.
- Chromosome sizes file.
- Bowtie2 index prefix.
- HiC-Pro restriction fragment BED.
- HiC-Pro 3 Singularity image.

Default HiC-Pro container:

```bash
/home/mowp/software/HiC-Pro/hicpro3.sif
```

## Quick start

```bash
workflows/hic_hichip/scripts/deploy_pipeline.sh /path/to/project_hichip
cd /path/to/project_hichip
# edit config/config.yml, config/samples.tsv, and config/comparisons.tsv
./run_snakemake.sh -n --cores 1
./run_snakemake.sh --cores 32 --rerun-incomplete --printshellcmds
```

## Configuration

Key sections in `config/config.yml`:

- `genome`: reference files.
- `preprocessing`: optional fastp settings. `enable_fastp` defaults to `false` for Hi-C/HiChIP.
- `hicpro`: Singularity image, HiC-Pro config template, resolutions, enzyme/ligation settings.
- `loop_calling`: FitHiChIP/Mustache/cooltools settings.
- `diffloop`: loop universe, quantification source, DESeq2 thresholds.
- `qc`: MultiQC and contact QC options.
- `annotation`: optional promoter/enhancer/GTF resources.

## Output structure

```text
results/
├── clean_data/
├── hicpro_input/
├── hicpro/
├── valid_pairs/
├── matrix/
├── cool/
├── qc/
├── loops/
│   └── consensus/
├── diffloop/
│   └── counts/
├── annotation/
├── reports/
├── logs/
└── benchmarks/
```

Important final outputs:

- `results/valid_pairs/<sample>.allValidPairs`
- `results/cool/<sample>.mcool`
- `results/loops/<sample>/<sample>.loops.bedpe`
- `results/loops/consensus/loop_universe.bedpe`
- `results/diffloop/counts/loop_counts.tsv`
- `results/diffloop/<comparison>/diffloops.tsv`
- `results/reports/hic_hichip_report.html`

## Loop calling strategy

For **HiChIP**, the recommended caller is `fithichip`. For biologically meaningful HiChIP loop calling, provide peak files, especially for marks such as H3K27ac, CTCF, RAD21, or RNAPII.

For **Hi-C**, `mustache` or `cooltools` is easier to deploy than GPU-dependent HiCCUPS. HiCCUPS can be added later as an optional rule if the compute environment supports it.

FitHiChIP is an external script package with Python/R/command-line dependencies. For a project-local install, run:

```bash
./scripts/install_fithichip.sh --prefix tools/FitHiChIP
```

Then set these values in `config/run_snakemake.env`:

```bash
FITHICHIP_HOME=tools/FitHiChIP
FITHICHIP_SCRIPT=tools/FitHiChIP/FitHiChIP_HiCPro.sh
```

You can also import precomputed loop files by adding this to `config/config.yml`:

```yaml
loop_calling:
  precomputed_loops:
    sample1: /path/to/sample1.loops.bedpe
    sample2: /path/to/sample2.loops.bedpe
```

## Differential loop interpretation

The default differential loop table uses this direction:

- `log2FoldChange > 0`: group1 has stronger loop contact than group2.
- `log2FoldChange < 0`: group2 has stronger loop contact than group1.

Formal DESeq2 statistics require biological replicates. Without sufficient replicates, the workflow writes descriptive fold-change results and p-values/FDR are set to `NA`.

## Biological notes

Do not interpret all loops as enhancer-promoter loops.

- H3K27ac HiChIP is enriched for active regulatory contacts and often helps prioritize enhancer/promoter loops.
- CTCF or RAD21 HiChIP is more architectural/cohesin-related.
- Hi-C loops are not mark-specific and require separate annotation to infer regulatory function.

## Troubleshooting

### HiC-Pro container not found

Check:

```yaml
hicpro:
  container: /home/mowp/software/HiC-Pro/hicpro3.sif
```

### Missing HiC-Pro matrix or validPairs

Check `results/logs/hicpro/run_hicpro.log` and confirm that HiC-Pro completed successfully. Also confirm that FASTQ filenames in `results/hicpro_input/<sample>/` match `PAIR1_EXT` and `PAIR2_EXT` in the generated config.

### FitHiChIP fails

Run `./scripts/install_fithichip.sh --prefix tools/FitHiChIP`, then check that `FITHICHIP_SCRIPT` in `config/run_snakemake.env` or `loop_calling.fithichip_script` in `config/config.yml` points to the runner. FitHiChIP also requires peaks through `samples.tsv` or `loop_calling.external_peaks`.

### Differential loops have NA p-values

This usually means there are not enough biological replicates or DESeq2 is unavailable. The output can still be used descriptively, but not as a formal statistical test.

## Notes for AI agents

Before running, inspect:

- `README.md`
- `agent/manifest.yml`
- `agent/io.contract.yml`
- `config/config.yml`
- `config/samples.tsv`
- `config/comparisons.tsv`

Always run a dry-run before a real execution:

```bash
./run_snakemake.sh -n --cores 1
```
