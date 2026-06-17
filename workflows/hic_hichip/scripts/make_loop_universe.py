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
    p.add_argument("--group-consensus-min-replicates", type=int, default=1)
    p.add_argument("--group-consensus-min-fraction", type=float, default=0.5)
    p.add_argument("--universe", required=True)
    p.add_argument("--group-consensus", required=True)
    p.add_argument("--annotation", required=True)
    return p.parse_args()


def read_loop_table(path):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=BEDPE_COLUMNS)
    missing = [c for c in BEDPE_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"Loop table missing columns {missing}: {path}")
    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=BEDPE_COLUMNS)


def canonical_anchor_pair(a, b):
    return (a, b) if a <= b else (b, a)


def cluster_anchors(anchors, slop):
    if not anchors:
        return {}, {}
    df = pd.DataFrame(anchors, columns=["anchor_key", "chrom", "start", "end"]).drop_duplicates()
    df["start_s"] = (df["start"].astype(int) - slop).clip(lower=0)
    df["end_s"] = df["end"].astype(int) + slop
    mapping = {}
    records = {}
    cluster_id = 0
    for chrom, sub in df.sort_values(["chrom", "start_s", "end_s"]).groupby("chrom", sort=False):
        current_start = None
        current_end = None
        current_keys = []
        for row in sub.itertuples(index=False):
            if current_end is None or int(row.start_s) > current_end:
                if current_keys:
                    cid = f"anchor_{cluster_id}"
                    starts = [k[1] for k in current_keys]
                    ends = [k[2] for k in current_keys]
                    records[cid] = (chrom, min(starts), max(ends))
                    for key in current_keys:
                        mapping[key] = cid
                    cluster_id += 1
                current_start = int(row.start_s)
                current_end = int(row.end_s)
                current_keys = [row.anchor_key]
            else:
                current_end = max(current_end, int(row.end_s))
                current_keys.append(row.anchor_key)
        if current_keys:
            cid = f"anchor_{cluster_id}"
            starts = [k[1] for k in current_keys]
            ends = [k[2] for k in current_keys]
            records[cid] = (chrom, min(starts), max(ends))
            for key in current_keys:
                mapping[key] = cid
            cluster_id += 1
    return mapping, records


def support_summary(values):
    return ",".join(sorted(set(map(str, values))))


def group_support(grouped, all_groups, group_sizes):
    support = {}
    for group in all_groups:
        n = len(set(grouped.loc[grouped["group"] == group, "sample"]))
        support[group] = n
    return ";".join(f"{g}:{support[g]}/{group_sizes.get(g, 0)}" for g in sorted(support))


def main():
    args = parse_args()
    sample_names = [s for s in args.samples.split(",") if s]
    sample_table = pd.read_csv(args.sample_table, sep="\t")
    group_of = dict(zip(sample_table["sample"], sample_table["group"]))
    group_sizes = sample_table.groupby("group")["sample"].nunique().to_dict()
    all_groups = sorted(group_sizes)

    loop_rows = []
    anchors = []
    for sample, path in zip(sample_names, args.loops):
        df = read_loop_table(path)
        for row in df.itertuples(index=False):
            a = (str(row.chrom1), int(row.start1), int(row.end1))
            b = (str(row.chrom2), int(row.start2), int(row.end2))
            a, b = canonical_anchor_pair(a, b)
            loop_rows.append({"sample": sample, "group": group_of.get(sample, "NA"), "a": a, "b": b})
            anchors.append((a, *a))
            anchors.append((b, *b))

    anchor_to_cluster, cluster_records = cluster_anchors(anchors, args.anchor_slop)
    clustered_rows = []
    for row in loop_rows:
        c1, c2 = canonical_anchor_pair(anchor_to_cluster[row["a"]], anchor_to_cluster[row["b"]])
        clustered_rows.append({**row, "anchor1_cluster": c1, "anchor2_cluster": c2})

    if clustered_rows:
        all_loops = pd.DataFrame(clustered_rows).drop_duplicates(["sample", "anchor1_cluster", "anchor2_cluster"])
    else:
        all_loops = pd.DataFrame(columns=["sample", "group", "anchor1_cluster", "anchor2_cluster"])

    records = []
    for (c1, c2), sub in all_loops.groupby(["anchor1_cluster", "anchor2_cluster"], dropna=False):
        a = cluster_records[c1]
        b = cluster_records[c2]
        sample_support = set(sub["sample"])
        groups = set(sub["group"])
        support_by_group = group_support(sub, all_groups, group_sizes)
        passes_group = False
        for group in all_groups:
            n = len(set(sub.loc[sub["group"] == group, "sample"]))
            required = max(args.group_consensus_min_replicates, int(args.group_consensus_min_fraction * group_sizes.get(group, 0) + 0.999999))
            if n >= required:
                passes_group = True
        records.append({
            "chrom1": a[0], "start1": a[1], "end1": a[2],
            "chrom2": b[0], "start2": b[1], "end2": b[2],
            "n_samples": len(sample_support),
            "samples": support_summary(sample_support),
            "n_groups": len(groups),
            "groups": support_summary(groups),
            "support_by_group": support_by_group,
            "anchor1_cluster": c1,
            "anchor2_cluster": c2,
            "passes_group_consensus": passes_group,
        })
    counts = pd.DataFrame(records)
    if counts.empty:
        counts = pd.DataFrame(columns=BEDPE_COLUMNS + ["n_samples", "samples", "n_groups", "groups", "support_by_group", "anchor1_cluster", "anchor2_cluster", "passes_group_consensus"])

    if args.mode == "intersection":
        universe = counts[counts["n_samples"] == len(sample_names)].copy()
    elif args.mode == "group_consensus":
        universe = counts[counts["passes_group_consensus"]].copy()
    else:
        universe = counts.copy()

    universe = universe.sort_values(BEDPE_COLUMNS).reset_index(drop=True)
    universe["loop_id"] = [f"loop_{i+1}" for i in range(len(universe))]
    out_cols = BEDPE_COLUMNS + ["loop_id", "n_samples", "samples", "n_groups", "groups", "support_by_group", "anchor1_cluster", "anchor2_cluster"]
    for path in [args.universe, args.group_consensus, args.annotation]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    universe[out_cols].to_csv(args.universe, sep="\t", index=False)
    counts[counts["passes_group_consensus"]].assign(
        loop_id=lambda x: [f"group_loop_{i+1}" for i in range(len(x))]
    )[out_cols].to_csv(args.group_consensus, sep="\t", index=False)
    annot = universe[out_cols].copy()
    annot["loop_class"] = annot.apply(lambda r: "interchromosomal" if r.chrom1 != r.chrom2 else "intrachromosomal", axis=1)
    annot["distance"] = annot.apply(lambda r: abs(((int(r.start1)+int(r.end1))//2) - ((int(r.start2)+int(r.end2))//2)) if r.chrom1 == r.chrom2 else "NA", axis=1)
    annot.to_csv(args.annotation, sep="\t", index=False)


if __name__ == "__main__":
    main()
