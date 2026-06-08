#!/usr/bin/env python3
"""Convert HiC-Pro matrix/bin files to a cooler file."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--bins", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--resolution", required=True, type=int)
    p.add_argument("--output", required=True)
    p.add_argument("--balance", default="true")
    return p.parse_args()


def read_bins(path: str) -> pd.DataFrame:
    bins = pd.read_csv(path, sep="\t", header=None, comment="#")
    # HiC-Pro abs.bed is typically: chrom start end abs_id
    if bins.shape[1] < 4:
        raise ValueError(f"Expected at least 4 columns in HiC-Pro abs bed: {path}")
    bins = bins.iloc[:, :4]
    bins.columns = ["chrom", "start", "end", "abs_id"]
    bins["start"] = bins["start"].astype(int)
    bins["end"] = bins["end"].astype(int)
    bins["abs_id"] = bins["abs_id"].astype(int)
    bins = bins.sort_values("abs_id").reset_index(drop=True)
    bins["bin_id"] = range(len(bins))
    return bins


def read_pixels(path: str, bins: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", header=None, comment="#")
    if matrix.shape[1] < 3:
        raise ValueError(f"Expected 3 columns in HiC-Pro sparse matrix: {path}")
    matrix = matrix.iloc[:, :3]
    matrix.columns = ["bin1_abs", "bin2_abs", "count"]
    mapper = dict(zip(bins["abs_id"], bins["bin_id"]))
    matrix["bin1_id"] = matrix["bin1_abs"].map(mapper)
    matrix["bin2_id"] = matrix["bin2_abs"].map(mapper)
    matrix = matrix.dropna(subset=["bin1_id", "bin2_id"])
    matrix["bin1_id"] = matrix["bin1_id"].astype(int)
    matrix["bin2_id"] = matrix["bin2_id"].astype(int)
    matrix["count"] = pd.to_numeric(matrix["count"], errors="coerce").fillna(0)
    matrix = matrix[matrix["count"] > 0]
    return matrix[["bin1_id", "bin2_id", "count"]]


def main():
    args = parse_args()
    try:
        import cooler
    except ImportError as exc:
        raise SystemExit("The 'cooler' Python package is required for matrix conversion.") from exc
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    bins = read_bins(args.bins)
    pixels = read_pixels(args.matrix, bins)
    cooler_bins = bins[["chrom", "start", "end"]].copy()
    cooler.create_cooler(str(out), cooler_bins, pixels, ordered=True, dtypes={"count": "float64"})
    if str(args.balance).lower() in {"true", "1", "yes", "y"}:
        try:
            cooler.balance_cooler(cooler.Cooler(str(out)), store=True)
        except Exception as exc:  # keep the unbalanced cool usable
            print(f"WARNING: cooler balancing failed: {exc}")


if __name__ == "__main__":
    main()
