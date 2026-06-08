# Agent task templates for BioFlowAI hic_hichip

## Deploy this workflow
Read `workflows/hic_hichip/README.md` and `agent/manifest.yml`, then deploy the workflow into my project directory. Do not modify raw FASTQ or reference files. Prepare `config/config.yml`, `samples.tsv`, and `comparisons.tsv`, then run a dry-run only.

## Fill config
Update `config/config.yml` using my project paths. Confirm the HiC-Pro Singularity image is `/home/mowp/software/HiC-Pro/hicpro3.sif`. Check genome FASTA, chrom sizes, Bowtie2 index prefix, and restriction fragment BED.

## Dry-run
Run `./run_snakemake.sh dry-run 8`. If it fails, report the missing input or invalid config key and suggest the smallest fix.

## Debug failed HiC-Pro
Inspect `results/logs/hicpro/run_hicpro.log`, the generated HiC-Pro config, and the structure under `results/hicpro_input/`. Check whether pair suffixes, reference genome, restriction fragments, and bind paths are correct.

## Debug missing loops
Check `results/logs/loops/call_loops/<sample>.log`. For HiChIP FitHiChIP, confirm peak BED is provided. For Hi-C, consider setting `loop_calling.caller` to `mustache` or `cooltools`.

## Explain diffloop result
Read `results/diffloop/<comparison>/diffloops.tsv`. Explain that positive log2FoldChange means group1 has stronger loop contact than group2. Highlight significant loops by padj, effect size, and whether they connect promoter/enhancer annotations if available.
