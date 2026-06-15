#!/usr/bin/env python3
"""Format bedtools multicov output as a peak-by-sample count matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Format bedtools multicov output.")
    parser.add_argument("--multicov", required=True, help="Input bedtools multicov TSV.")
    parser.add_argument("--samples", nargs="+", required=True, help="Sample names in BAM order.")
    parser.add_argument("--output", required=True, help="Output count matrix TSV.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.multicov).open(newline="") as handle, Path(args.output).open("w", newline="") as output:
        reader = csv.reader(handle, delimiter="\t")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["peak_id", "chrom", "start", "end", *args.samples])
        for row in reader:
            if len(row) < 6 + len(args.samples):
                raise ValueError("multicov row has fewer columns than expected.")
            writer.writerow([row[3], row[0], row[1], row[2], *row[6 : 6 + len(args.samples)]])


if __name__ == "__main__":
    main()
