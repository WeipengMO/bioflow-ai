#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--loops", required=True)
    p.add_argument("--promoters", default="")
    p.add_argument("--enhancers", default="")
    p.add_argument("--gtf", default="")
    p.add_argument("--loops-to-genes", required=True)
    p.add_argument("--promoter-enhancer", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    loops = pd.read_csv(args.loops, sep="\t")
    # Lightweight placeholder annotation: preserve loop table and flag whether resources were supplied.
    # For production, replace this with bedtools intersect against promoter/enhancer BED and GTF-derived promoters.
    ann = loops.copy()
    ann["promoter_resource"] = args.promoters or args.gtf or "not_provided"
    ann["enhancer_resource"] = args.enhancers or "not_provided"
    ann["nearest_gene"] = "NA"
    ann["annotation_note"] = "Run bedtools-based promoter/enhancer annotation by supplying annotation.promoters/enhancers/gtf."
    pe = ann.copy()
    pe["loop_type"] = "unclassified"
    for p in [args.loops_to_genes, args.promoter_enhancer]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    ann.to_csv(args.loops_to_genes, sep="\t", index=False)
    pe.to_csv(args.promoter_enhancer, sep="\t", index=False)


if __name__ == "__main__":
    main()
