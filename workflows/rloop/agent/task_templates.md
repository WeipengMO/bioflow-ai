# Agent Task Templates

## Deploy a project directory

```bash
cd workflows/rloop
./scripts/deploy_pipeline.sh /path/to/rloop_project
```

Use `--force` to replace deploy-managed files and links when re-deploying into an existing directory:

```bash
./scripts/deploy_pipeline.sh /path/to/rloop_project --force
```

After deployment:

1. Put FASTQ files under `raw_data/`.
2. Edit `config/config.yml` and `config/samples.yml`.
3. Dry-run with `./run_snakemake.sh -n --cores 1`.

## Dry-run the workflow

```bash
cd workflows/rloop
snakemake -np --use-conda
```

## Configure a DRIP-seq project with input and RNase H controls

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

Set broad peak mode unless the experiment expects sharp sites:

```yaml
peak_type: broad
enable_rnaseh_subtraction: true  # writes rnaseh_no_overlap plus deprecated rnaseh_sensitive alias
```

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

## Set up spike-in normalization

1. Run `bash scripts/build_spikein_index.sh` to download E. coli K-12 MG1655 and build a Bowtie2 index.
2. In `config/config.yml`, set:
   ```yaml
   signal_scale_factor_method: spikein
   spikein_genome: data/spikein/ecoli_mg1655
   spikein_method: ratio
   spikein_reference_sample: <reference_sample_name>
   enable_common_scale_bigwig: true
   ```
3. Dry-run with `snakemake -np --use-conda`.

## Set up CUT&Tag duplicate handling

Default `peak_duplicate_mode: auto` uses markdup BAM for peak calling with `--keep-dup all`.
Only switch to `peak_duplicate_mode: remove` when Picard duplication rate is extremely high (80-90%+):

```yaml
peak_duplicate_mode: auto  # auto | keep_marked | remove
dup_threshold_for_dedup: 0.90
```
