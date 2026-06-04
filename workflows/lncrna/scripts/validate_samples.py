#!/usr/bin/env python3
"""Validate sample sheet for lncrna_discovery workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("")
    required = {"sample", "fq1", "fq2"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")
    if df["sample"].duplicated().any():
        raise SystemExit("Duplicated sample names: " + ", ".join(df.loc[df["sample"].duplicated(), "sample"]))
    missing_files = []
    for _, row in df.iterrows():
        for col in ["fq1", "fq2"]:
            path = row[col]
            if path and not Path(path).exists():
                missing_files.append(f"{row['sample']}:{col}:{path}")
    if missing_files:
        raise SystemExit("Missing FASTQ files:\n" + "\n".join(missing_files))
    print(f"Sample sheet OK: {len(df)} samples")


if __name__ == "__main__":
    main()
