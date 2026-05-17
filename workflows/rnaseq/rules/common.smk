import os
import pandas as pd
from pathlib import Path

OUTDIR = config.get("outdir", "results")
SAMPLES_TSV = config["samples"]
ENV = "../envs/rnaseq_expression.yml"

samples_df = pd.read_csv(SAMPLES_TSV, sep="\t", dtype=str).fillna("")
required_cols = {"sample", "fq1", "fq2"}
missing_cols = required_cols - set(samples_df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns in {SAMPLES_TSV}: {sorted(missing_cols)}")

samples_df["sample"] = samples_df["sample"].astype(str)
if samples_df["sample"].duplicated().any():
    duplicated = samples_df.loc[samples_df["sample"].duplicated(), "sample"].tolist()
    raise ValueError(f"Duplicated sample names: {duplicated}")

SAMPLES = samples_df["sample"].tolist()
SAMPLE_TO_FQ1 = dict(zip(samples_df["sample"], samples_df["fq1"]))
SAMPLE_TO_FQ2 = dict(zip(samples_df["sample"], samples_df["fq2"]))

PE_SAMPLES = [s for s in SAMPLES if SAMPLE_TO_FQ2.get(s, "") not in ["", "nan", "None"]]
SE_SAMPLES = [s for s in SAMPLES if s not in PE_SAMPLES]


def is_paired(sample):
    return sample in PE_SAMPLES


def raw_fq1(wildcards):
    return SAMPLE_TO_FQ1[wildcards.sample]


def raw_fq2(wildcards):
    return SAMPLE_TO_FQ2[wildcards.sample]


def clean_reads(wildcards):
    sample = wildcards.sample
    if is_paired(sample):
        return [
            f"{OUTDIR}/clean_fastq/{sample}.R1.clean.fastq.gz",
            f"{OUTDIR}/clean_fastq/{sample}.R2.clean.fastq.gz",
        ]
    return [f"{OUTDIR}/clean_fastq/{sample}.clean.fastq.gz"]


def featurecounts_extra(wildcards):
    fc = config.get("featurecounts", {})
    base = fc.get("extra", "")
    layout_extra = fc.get("paired_extra", "-p -B -C") if is_paired(wildcards.sample) else fc.get("single_extra", "")
    return " ".join(x for x in [layout_extra, base] if x)


def strand_flag_featurecounts():
    strandness = config.get("strandness", "unstranded")
    mapping = {"unstranded": 0, "forward": 1, "reverse": 2}
    if strandness not in mapping:
        raise ValueError(f"Invalid strandness: {strandness}. Use unstranded, forward, or reverse.")
    return mapping[strandness]
