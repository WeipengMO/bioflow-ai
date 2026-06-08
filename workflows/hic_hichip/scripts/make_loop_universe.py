#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

BEDPE_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--loops", nargs="+", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--sample-table", required=True)
    p.add_argument("--mode", default="union", choices=["union", "intersection", "group_consensus"])
    p.add_argument("--anchor-slop", type=int, default=0)
    p.add_argument("--universe", required=True)
    p.add_argument("--group-consensus", required=True)
    p.add_argument("--annotation", required=True)
    return p.parse_args()


def canonicalize(row, slop=0):
    a = (str(row.chrom1), max(0, int(row.start1)-slop), int(row.end1)+slop)
    b = (str(row.chrom2), max(0, int(row.start2)-slop), int(row.end2)+slop)
    if a > b:
        a, b = b, a
    return (*a, *b)


def read_sample_table(path):
    return pd.read_csv(path, sep="\t")


def main():
    args = parse_args()
    sample_names = [s for s in args.samples.split(",") if s]
    sample_table = read_sample_table(args.sample_table)
    group_of = dict(zip(sample_table["sample"], sample_table["group"]))
    frames = []
    for sample, path in zip(sample_names, args.loops):
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        df = df.dropna(subset=BEDPE_COLUMNS)
        keys = df.apply(lambda r: canonicalize(r, args.anchor_slop), axis=1)
        norm = pd.DataFrame(list(keys), columns=BEDPE_COLUMNS)
        norm["sample"] = sample
        norm["group"] = group_of.get(sample, "NA")
        frames.append(norm)
    if frames:
        all_loops = pd.concat(frames, ignore_index=True).drop_duplicates()
    else:
        all_loops = pd.DataFrame(columns=BEDPE_COLUMNS + ["sample", "group"])
    counts = all_loops.groupby(BEDPE_COLUMNS, dropna=False).agg(
        n_samples=("sample", "nunique"),
        samples=("sample", lambda x: ",".join(sorted(set(map(str, x))))),
        n_groups=("group", "nunique"),
        groups=("group", lambda x: ",".join(sorted(set(map(str, x)))))
    ).reset_index()
    if args.mode == "intersection":
        universe = counts[counts["n_samples"] == len(sample_names)].copy()
    elif args.mode == "group_consensus":
        universe = counts[counts["n_groups"] >= 1].copy()
    else:
        universe = counts.copy()
    universe = universe.sort_values(BEDPE_COLUMNS).reset_index(drop=True)
    universe["loop_id"] = [f"loop_{i+1}" for i in range(len(universe))]
    out_cols = BEDPE_COLUMNS + ["loop_id", "n_samples", "samples", "n_groups", "groups"]
    for path in [args.universe, args.group_consensus, args.annotation]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    universe[out_cols].to_csv(args.universe, sep="\t", index=False)
    # group consensus: keep loops observed in at least one sample per present group; this is an annotation table for now.
    universe[out_cols].to_csv(args.group_consensus, sep="\t", index=False)
    annot = universe[out_cols].copy()
    annot["loop_class"] = annot.apply(lambda r: "interchromosomal" if r.chrom1 != r.chrom2 else "intrachromosomal", axis=1)
    annot["distance"] = annot.apply(lambda r: abs(((int(r.start1)+int(r.end1))//2) - ((int(r.start2)+int(r.end2))//2)) if r.chrom1 == r.chrom2 else "NA", axis=1)
    annot.to_csv(args.annotation, sep="\t", index=False)


if __name__ == "__main__":
    main()
