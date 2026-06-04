# lncRNA Discovery Workflow

This workflow follows the repository's RNA-seq workflow layout and reuses the same preprocessing/alignment/counting environment wherever the tools are identical. It adds transcript assembly, annotation comparison, candidate lncRNA filtering, FASTA extraction, and lncRNA expression matrix generation.

## Inputs

Project files:

```text
config/config.yml
config/samples.tsv
raw_data/*.fastq.gz or *.fq.gz
```

`samples.tsv` example:

```tsv
sample	fq1	fq2	group	replicate
control_1	raw_data/control_1_R1.fastq.gz	raw_data/control_1_R2.fastq.gz	control	1
control_2	raw_data/control_2_R1.fastq.gz	raw_data/control_2_R2.fastq.gz	control	2
treated_1	raw_data/treated_1_R1.fastq.gz	raw_data/treated_1_R2.fastq.gz	treated	1
treated_2	raw_data/treated_2_R1.fastq.gz	raw_data/treated_2_R2.fastq.gz	treated	2
```

- `sample`: unique sample name (required)
- `fq1`: R1 FASTQ or single-end FASTQ (required)
- `fq2`: R2 FASTQ for paired-end data; leave empty for single-end data (optional)
- `group`: sample group for downstream differential analysis (optional)
- `replicate`: replicate ID for distinguishing biological replicates (optional)
- Extra columns are allowed for downstream use

### Reference Resource Preparation

All genome-related files must use the same assembly.

```yaml
genome:
  fasta: path/to/hg38.fa
  gtf: path/to/hg38.gtf
  star_index: path/to/hg38_star_index
```

Build the STAR index:

```bash
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir path/to/hg38_star_index \
  --genomeFastaFiles path/to/hg38.fa \
  --sjdbGTFfile path/to/hg38.gtf
```

## Configuration

Main config options:

```yaml
genome:
  fasta: /path/to/genome.fa
  gtf: /path/to/annotation.gtf
  star_index: /path/to/star_index

strandness: unstranded  # unstranded / forward / reverse

lncrna_filter:
  min_length: 200
  min_exons: 2
  class_codes: [u, i, x, o, e]

coding_filter:
  max_orf_aa: 100
  exclude_transcripts: ""
```

- `strandness` controls featureCounts `-s` parameter:
  - `unstranded` -> `-s 0`
  - `forward` -> `-s 1`
  - `reverse` -> `-s 2`
- `lncrna_filter.class_codes` uses `gffcompare` classes:
  - `u`: intergenic
  - `i`: intronic
  - `x`: antisense
  - `o`: generic exonic overlap
  - `e`: single-exon transfrag overlap

## Workflow Logic

1. `fastp`: raw FASTQ QC and filtering; uses the same conda environment as RNA-seq.
2. `STAR`: splice-aware alignment to sorted BAM; uses the same conda environment as RNA-seq.
3. `StringTie`: per-sample reference-guided transcript assembly.
4. `StringTie --merge`: cross-sample transcript merge.
5. `gffcompare`: compare merged transcripts against the reference annotation and assign class codes.
6. `filter_lncRNA_candidates.py`: retain candidate lncRNAs by class code, exon-union length, and exon count.
7. `gffread`: extract candidate lncRNA transcript sequences.
8. `filter_high_confidence_lncRNA.py`: remove candidates with long six-frame ORFs and optionally remove transcript IDs flagged by external coding-potential tools.
9. `featureCounts`: quantify reads over both the reference gene annotation and high-confidence lncRNA exons; uses the same conda environment as RNA-seq.
10. `matrix`: merge raw counts and compute whole-transcriptome gene matrices plus high-confidence lncRNA matrices.
11. `bamCoverage` and `MultiQC`: generate bigWig tracks and QC summary; uses the same conda environment as RNA-seq.

## Coding Potential Filtering

By default, the workflow applies an ORF length filter:

```yaml
coding_filter:
  max_orf_aa: 100
```

Transcripts whose longest six-frame ORF is greater than `max_orf_aa` are removed from the high-confidence lncRNA set. If CPAT, CPC2, PLEK, Pfam/HMMER, or another external method has already produced a newline-delimited list of transcript IDs predicted to be protein-coding, provide it with:

```yaml
coding_filter:
  exclude_transcripts: path/to/coding_transcript_ids.txt
```

The final quantification uses `high_confidence_lncRNA.gtf`, not the broader `candidate_lncRNA.gtf`.

## Environment Reuse

The reusable RNA-seq environment file is intentionally kept identical to:

```text
workflows/rnaseq/envs/rnaseq_expression.yml
```

This lets Snakemake reuse the same conda environment for shared tools such as `fastp`, `STAR`, `samtools`, `featureCounts`, `deepTools`, and `MultiQC`. lncRNA-specific discovery tools are isolated in:

```text
workflows/lncrna/envs/lncrna_discovery.yml
```

## Outputs

```text
results/lncrna/candidate_lncRNA.gtf          # Filtered candidate lncRNA annotation
results/lncrna/candidate_lncRNA.fa           # Candidate lncRNA transcript FASTA
results/lncrna/candidate_lncRNA.summary.tsv  # Candidate class, length, exon count summary
results/lncrna/high_confidence_lncRNA.gtf    # Candidate set after coding-potential filtering
results/lncrna/high_confidence_lncRNA.fa
results/lncrna/high_confidence_lncRNA.summary.tsv
results/lncrna/coding_filter/orf_metrics.tsv # Longest ORF metrics and keep/drop calls
results/matrix/gene_counts.tsv               # Whole-transcriptome/gene raw count matrix
results/matrix/gene_fpkm.tsv                 # Whole-transcriptome/gene FPKM matrix
results/matrix/gene_tpm.tsv                  # Whole-transcriptome/gene TPM matrix
results/matrix/gene_length.tsv
results/matrix/lncrna_counts.tsv             # High-confidence lncRNA raw count matrix
results/matrix/lncrna_fpkm.tsv               # High-confidence lncRNA FPKM matrix
results/matrix/lncrna_tpm.tsv                # High-confidence lncRNA TPM matrix
results/matrix/lncrna_length.tsv
results/bam/{sample}.sorted.bam              # Aligned BAM
results/bam/{sample}.sorted.bam.bai
results/bigwig/{sample}.bw                   # bigWig tracks
results/gffcompare/merged.stats              # Annotation comparison statistics
results/report/multiqc_report.html           # QC summary
```

## Quick Start

Create a project directory from this workflow template:

```bash
cd workflows/lncrna
./scripts/deploy_pipeline.sh /path/to/lncrna_project
cd /path/to/lncrna_project
```

Place FASTQ files under `raw_data/`, edit `config/config.yml` and `config/samples.tsv`, then run:

```bash
./run_snakemake.sh -n
./run_snakemake.sh
```

## Notes for AI Agents

- Read this README and the files under `agent/` before modifying config or starting a run.
- Prefer deploying the workflow into a separate project directory instead of editing files under `workflows/lncrna/` directly.
- Preserve existing `raw_data/` files and treat user-provided sample tables and config files as user-owned inputs unless the user explicitly requests edits.
- Start with `./run_snakemake.sh -n -p` before attempting a full run.
- When reporting back, summarize candidate annotation files, expression matrices, unresolved required inputs, and the next recommended command.

## External Resource Links

- STAR: https://github.com/alexdobin/STAR
- StringTie: https://ccb.jhu.edu/software/stringtie/
- gffcompare/gffread: https://ccb.jhu.edu/software/stringtie/gff.shtml
- Subread/featureCounts: http://subread.sourceforge.net/
- deepTools: https://deeptools.readthedocs.io/
- MultiQC: https://multiqc.info/

For optional downstream analysis, use DESeq2 or edgeR for differential lncRNA expression and combine candidates with coding-potential tools such as CPAT/CPC2 if a stricter novel-lncRNA discovery study is required.
