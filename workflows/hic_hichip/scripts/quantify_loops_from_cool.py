#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--universe", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--sample-table", required=True)
    p.add_argument("--cools", nargs="+", required=True)
    p.add_argument("--resolution", type=int, required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--loop-metadata", required=True)
    p.add_argument("--sample-metadata", required=True)
    p.add_argument("--threads", type=int, default=1)  # reserved
    p.add_argument("--balanced", default="false")
    return p.parse_args()


def region(chrom, start, end):
    return f"{chrom}:{int(start)}-{int(end)}"


def as_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes", "y", "on"}


def loop_count_from_cool(cool_path, loops, balanced=False):
    try:
        import cooler
    except ImportError as exc:
        raise SystemExit("The cooler Python package is required for cool-based loop quantification.") from exc
    c = cooler.Cooler(cool_path)
    if loops.empty:
        return []
    out = []
    available_chroms = set(c.chromnames)
    matrix = c.matrix(balance=balanced)
    for _, row in loops.iterrows():
        try:
            if row.chrom1 not in available_chroms or row.chrom2 not in available_chroms:
                print(f"WARNING: loop chrom missing from cooler {cool_path}: {row.chrom1}, {row.chrom2}", file=sys.stderr)
                out.append(0.0)
                continue
            sub = matrix.fetch(region(row.chrom1, row.start1, row.end1), region(row.chrom2, row.start2, row.end2))
            out.append(float(sub.sum()))
        except Exception as exc:
            print(f"WARNING: failed to quantify loop from {cool_path}: {exc}", file=sys.stderr)
            out.append(0.0)
    return out


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    loops = pd.read_csv(args.universe, sep="\t")
    if "loop_id" not in loops.columns:
        loops["loop_id"] = [f"loop_{i+1}" for i in range(len(loops))]
    counts = pd.DataFrame({"loop_id": loops["loop_id"]})
    for sample, cool in zip(samples, args.cools):
        counts[sample] = loop_count_from_cool(cool, loops, balanced=as_bool(args.balanced))
    sample_meta = pd.read_csv(args.sample_table, sep="\t")
    Path(args.counts).parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.counts, sep="\t", index=False)
    loops.to_csv(args.loop_metadata, sep="\t", index=False)
    sample_meta.to_csv(args.sample_metadata, sep="\t", index=False)


if __name__ == "__main__":
    main()
