# Agent Task Templates

## Deploy a project directory

```bash
cd workflows/cut_tag_rloop
./scripts/deploy_pipeline.sh /path/to/cut_tag_rloop_project
```

Use `--force` to replace deploy-managed files and links when re-deploying into an existing directory:

```bash
./scripts/deploy_pipeline.sh /path/to/cut_tag_rloop_project --force
```

After deployment:

1. Put FASTQ files under `raw_data/`.
2. Edit `config/config.yml` and `config/samples.yml`.
3. Dry-run with `./run_snakemake.sh -n --cores 1`.

## Dry-run the workflow

```bash
cd workflows/cut_tag_rloop
snakemake -np --use-conda
```

## Configure CUT&Tag R-loop samples

```yaml
WT_CutTag_RLoop_rep1:
  role: treatment
  group: WT
  control: WT_Input_rep1
  rnaseh_control: WT_RNaseH_rep1
WT_Input_rep1:
  role: control
WT_RNaseH_rep1:
  role: rnaseh_control
```

Recommended defaults:

```yaml
assay: CUT&Tag-R-loop
peak_type: broad
scale_methods:
  - cpm
  - absolute_spikein
  - matched_ref_spikein
spikein:
  enabled: true
  genome: ecoli
  bowtie2_index: data/spikein/ecoli_mg1655
rnaseh_sensitive:
  mode: both
  min_fold_change: 2.0
```

## Confirm default spike-in scale

1. Run `bash scripts/build_spikein_index.sh` when `data/spikein/ecoli_mg1655` is not present.
2. In `config/config.yml`, confirm:
   ```yaml
   scale_methods:
     - cpm
     - absolute_spikein
     - matched_ref_spikein
   spikein:
     enabled: true
     genome: ecoli
     bowtie2_index: data/spikein/ecoli_mg1655
     min_spikein_reads: 1000
     warn_low_fraction: 0.001
   ```
3. Dry-run with `snakemake -np --use-conda`.
4. After execution, inspect `results/qc/normalization/spikein_summary.tsv` and `warnings/cut_tag_rloop_spikein.warning.tsv`.

## Inspect failed rule

1. Read the Snakemake error.
2. Identify the rule name and sample wildcard.
3. Open the matching file under `results/logs/` by default, or under `<outdir>/logs/` if the workflow config overrides `outdir`.
4. Check whether the failure came from missing inputs, inconsistent config, tool failure, or empty filtered output.

## Add a new sample

1. Add FASTQ files under `raw_data/` or update `config/fastq_manifest.tsv`.
2. Add the sample to `config/samples.yml`.
3. Set `role`, and for treatment samples set `group`, `control`, and `rnaseh_control` when applicable.
4. Run `snakemake -np --use-conda`.

## Debug FASTQ pairing

```bash
python scripts/make_fastq_manifest.py --raw-data-dir raw_data --output config/fastq_manifest.tsv
```

Then inspect `config/fastq_manifest.tsv` and set:

```yaml
fastq_manifest: config/fastq_manifest.tsv
```

## CUT&Tag duplicate handling

Default `peak_duplicate_mode: auto` uses duplicate-marked BAMs for peak calling with `--keep-dup all`.
Only switch to `peak_duplicate_mode: remove` as a high-duplication sensitivity analysis:

```yaml
peak_duplicate_mode: auto  # auto | keep_marked | remove
```
