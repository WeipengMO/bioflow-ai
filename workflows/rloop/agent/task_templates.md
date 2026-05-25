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
enable_rnaseh_subtraction: true
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
