#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs and compute FPKM/TPM matrices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd


def read_featurecounts(path: str) -> pd.Series:
    sample = Path(path).name.replace(".featureCounts.txt", "")
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.empty:
        raise ValueError(f"Empty featureCounts file: {path}")
    # featureCounts uses the BAM path as the last count column. This workflow counts one BAM per file.
    count_col = df.columns[-1]
    series = df.set_index("Geneid")[count_col].astype(float)
    series.name = sample
    return series


def merge_counts(files: Iterable[str], sample_order: List[str]) -> pd.DataFrame:
    series_list = [read_featurecounts(f) for f in files]
    counts = pd.concat(series_list, axis=1).fillna(0)
    available = [s for s in sample_order if s in counts.columns]
    extra = [s for s in counts.columns if s not in available]
    return counts[available + extra]


def compute_fpkm_tpm(counts: pd.DataFrame, gene_length: pd.DataFrame):
    lengths = gene_length.set_index("gene_id")["length"].astype(float)
    common = counts.index.intersection(lengths.index)
    if len(common) == 0:
        raise ValueError("No overlapping gene_id values between counts and gene_length table.")

    counts = counts.loc[common].astype(float)
    lengths_kb = lengths.loc[common] / 1000.0
    valid = lengths_kb > 0
    counts = counts.loc[valid]
    lengths_kb = lengths_kb.loc[valid]

    library_size_million = counts.sum(axis=0) / 1_000_000.0
    fpkm = counts.div(lengths_kb, axis=0).div(library_size_million.replace(0, np.nan), axis=1).fillna(0)

    rpk = counts.div(lengths_kb, axis=0)
    scaling = rpk.sum(axis=0) / 1_000_000.0
    tpm = rpk.div(scaling.replace(0, np.nan), axis=1).fillna(0)
    return counts.astype(int), fpkm, tpm


def add_gene_symbols(matrix: pd.DataFrame, gene_length: pd.DataFrame) -> pd.DataFrame:
    # gene_length.tsv has columns: gene_id, gene_name, chromosome, strand, length
    if "gene_name" in gene_length.columns:
        symbols = gene_length.set_index("gene_id")["gene_name"]
        matrix = matrix.copy()
        matrix.insert(0, "gene_symbol", symbols.reindex(matrix.index).fillna("").values)
    else:
        matrix = matrix.copy()
        matrix.insert(0, "gene_symbol", "")
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate raw count, FPKM, and TPM matrices from featureCounts outputs.")
    parser.add_argument("--featurecounts", nargs="+", required=True, help="Per-sample featureCounts output files.")
    parser.add_argument("--gene-length", required=True, help="Gene length TSV from build_gene_length_table.py.")
    parser.add_argument("--sample-sheet", required=True, help="Sample sheet with sample column.")
    parser.add_argument("--out-counts", required=True)
    parser.add_argument("--out-fpkm", required=True)
    parser.add_argument("--out-tpm", required=True)
    args = parser.parse_args()

    sample_df = pd.read_csv(args.sample_sheet, sep="\t", dtype=str).fillna("")
    sample_order = sample_df["sample"].tolist()

    counts = merge_counts(args.featurecounts, sample_order)
    gene_length = pd.read_csv(args.gene_length, sep="\t")
    counts, fpkm, tpm = compute_fpkm_tpm(counts, gene_length)

    counts = add_gene_symbols(counts, gene_length)
    fpkm = add_gene_symbols(fpkm, gene_length)
    tpm = add_gene_symbols(tpm, gene_length)

    for output in [args.out_counts, args.out_fpkm, args.out_tpm]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    counts.to_csv(args.out_counts, sep="\t", index_label="gene_id")
    fpkm.to_csv(args.out_fpkm, sep="\t", index_label="gene_id", float_format="%.6f")
    tpm.to_csv(args.out_tpm, sep="\t", index_label="gene_id", float_format="%.6f")


if __name__ == "__main__":
    main()
