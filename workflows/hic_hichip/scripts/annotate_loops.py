#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--loops", required=True)
    p.add_argument("--promoters", default="")
    p.add_argument("--enhancers", default="")
    p.add_argument("--gtf", default="")
    p.add_argument("--promoter-upstream", type=int, default=2000)
    p.add_argument("--promoter-downstream", type=int, default=2000)
    p.add_argument("--loops-to-genes", required=True)
    p.add_argument("--promoter-enhancer", required=True)
    return p.parse_args()


def attrs_get(attrs, key):
    m = re.search(rf'{key} "([^"]+)"', attrs)
    return m.group(1) if m else ""


def promoters_from_gtf(path, upstream, downstream):
    rows = []
    with open(path, errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end, strand, attrs = parts[0], int(parts[3]) - 1, int(parts[4]), parts[6], parts[8]
            biotype = attrs_get(attrs, "gene_type") or attrs_get(attrs, "gene_biotype")
            if biotype and biotype not in {"protein_coding", "lncRNA", "lincRNA"}:
                continue
            tss = start if strand != "-" else end
            rows.append({
                "chrom": chrom,
                "start": max(0, tss - upstream),
                "end": tss + downstream,
                "gene_id": attrs_get(attrs, "gene_id"),
                "gene_name": attrs_get(attrs, "gene_name") or attrs_get(attrs, "gene_id"),
                "strand": strand,
                "tss": tss,
            })
    return pd.DataFrame(rows)


def read_bed(path, kind):
    if not path or not Path(path).exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "id", "gene_id", "gene_name", "strand", "tss"])
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "id", "gene_id", "gene_name", "strand", "tss"])
    out = df.iloc[:, : min(df.shape[1], 6)].copy()
    cols = ["chrom", "start", "end", "id", "score", "strand"][: out.shape[1]]
    out.columns = cols
    if "id" not in out:
        out["id"] = [f"{kind}_{i+1}" for i in range(len(out))]
    out["gene_id"] = out["id"]
    out["gene_name"] = out["id"]
    out["tss"] = ((out["start"].astype(int) + out["end"].astype(int)) // 2)
    if "strand" not in out:
        out["strand"] = "."
    return out


def overlaps(anchor, features):
    if features.empty:
        return pd.DataFrame()
    chrom, start, end = anchor
    sub = features[(features["chrom"] == chrom) & (features["start"].astype(int) < end) & (features["end"].astype(int) > start)].copy()
    if not sub.empty and "tss" in sub:
        center = (start + end) // 2
        sub["distance_to_tss"] = (sub["tss"].astype(int) - center).abs()
    return sub


def main():
    args = parse_args()
    loops = pd.read_csv(args.loops, sep="\t")
    promoters = read_bed(args.promoters, "promoter")
    if promoters.empty and args.gtf:
        promoters = promoters_from_gtf(args.gtf, args.promoter_upstream, args.promoter_downstream)
    enhancers = read_bed(args.enhancers, "enhancer")
    ltog = []
    pe_rows = []
    for row in loops.itertuples(index=False):
        loop_id = getattr(row, "loop_id", "")
        anchors = {
            "anchor1": (row.chrom1, int(row.start1), int(row.end1)),
            "anchor2": (row.chrom2, int(row.start2), int(row.end2)),
        }
        prom_hits = {}
        enh_hits = {}
        for name, anchor in anchors.items():
            prom_hits[name] = overlaps(anchor, promoters)
            enh_hits[name] = overlaps(anchor, enhancers)
            for hit in prom_hits[name].itertuples(index=False):
                ltog.append({
                    "loop_id": loop_id, "anchor": name, "gene_id": hit.gene_id,
                    "gene_name": hit.gene_name, "distance_to_tss": getattr(hit, "distance_to_tss", "NA"),
                    "anchor_type": "promoter"
                })
        for p_anchor, e_anchor in [("anchor1", "anchor2"), ("anchor2", "anchor1")]:
            if prom_hits[p_anchor].empty or enh_hits[e_anchor].empty:
                continue
            for p_hit in prom_hits[p_anchor].itertuples(index=False):
                for e_hit in enh_hits[e_anchor].itertuples(index=False):
                    pe_rows.append({
                        "loop_id": loop_id,
                        "promoter_anchor": p_anchor,
                        "enhancer_anchor": e_anchor,
                        "gene_id": p_hit.gene_id,
                        "gene_name": p_hit.gene_name,
                        "enhancer_id": e_hit.id,
                        "loop_type": "promoter-enhancer_candidate",
                        "distance": abs(((int(row.start1)+int(row.end1))//2) - ((int(row.start2)+int(row.end2))//2)) if row.chrom1 == row.chrom2 else "NA",
                        "peak_source": getattr(row, "peak_source", "NA"),
                    })
    for p in [args.loops_to_genes, args.promoter_enhancer]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ltog, columns=["loop_id", "anchor", "gene_id", "gene_name", "distance_to_tss", "anchor_type"]).to_csv(args.loops_to_genes, sep="\t", index=False)
    pd.DataFrame(pe_rows, columns=["loop_id", "promoter_anchor", "enhancer_anchor", "gene_id", "gene_name", "enhancer_id", "loop_type", "distance", "peak_source"]).to_csv(args.promoter_enhancer, sep="\t", index=False)


if __name__ == "__main__":
    main()
