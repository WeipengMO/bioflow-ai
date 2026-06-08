#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import gzip
import bisect
import pandas as pd

BEDPE_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--universe", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--sample-table", required=True)
    p.add_argument("--validpairs", nargs="+", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--loop-metadata", required=True)
    p.add_argument("--sample-metadata", required=True)
    p.add_argument("--threads", type=int, default=1)  # reserved
    return p.parse_args()


def opener(path):
    return gzip.open if str(path).endswith(".gz") else open


def build_index(loops):
    by_chrom = {}
    for idx, row in loops.iterrows():
        for side in [1]:
            chrom = row["chrom1"]
            by_chrom.setdefault(chrom, []).append((int(row["start1"]), int(row["end1"]), idx))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def candidate_anchor_hits(index, chrom, pos):
    hits = []
    items = index.get(chrom, [])
    # simple scan around sorted starts; robust for moderate loop sets
    i = bisect.bisect_right(items, (pos, 10**18, 10**18))
    j = i - 1
    while j >= 0 and items[j][0] <= pos:
        start, end, idx = items[j]
        if end > pos:
            hits.append(idx)
        if start < pos - 1_000_000:
            break
        j -= 1
    return hits


def count_validpairs(path, loops):
    index = build_index(loops)
    counts = [0] * len(loops)
    with opener(path)(path, "rt", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            # HiC-Pro allValidPairs: readID chr1 pos1 strand1 chr2 pos2 strand2 ...
            chrom1, pos1, chrom2, pos2 = parts[1], int(parts[2]), parts[4], int(parts[5])
            for idx in candidate_anchor_hits(index, chrom1, pos1):
                row = loops.iloc[idx]
                ok = (
                    row.chrom2 == chrom2 and int(row.start2) <= pos2 < int(row.end2)
                ) or (
                    row.chrom1 == chrom2 and int(row.start1) <= pos2 < int(row.end1) and row.chrom2 == chrom1 and int(row.start2) <= pos1 < int(row.end2)
                )
                if ok:
                    counts[idx] += 1
    return counts


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    loops = pd.read_csv(args.universe, sep="\t")
    if "loop_id" not in loops.columns:
        loops["loop_id"] = [f"loop_{i+1}" for i in range(len(loops))]
    counts = pd.DataFrame({"loop_id": loops["loop_id"]})
    for sample, vp in zip(samples, args.validpairs):
        counts[sample] = count_validpairs(vp, loops)
    sample_meta = pd.read_csv(args.sample_table, sep="\t")
    Path(args.counts).parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.counts, sep="\t", index=False)
    loops.to_csv(args.loop_metadata, sep="\t", index=False)
    sample_meta.to_csv(args.sample_metadata, sep="\t", index=False)


if __name__ == "__main__":
    main()
