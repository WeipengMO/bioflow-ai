#!/usr/bin/env python3
"""Build exon-union gene length table from a GTF file."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_attrs(attr: str) -> Dict[str, str]:
    out = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        out[match.group(1)] = match.group(2)
    return out


def merge_intervals(intervals: List[Tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return sum(end - start + 1 for start, end in merged)


def build_gene_lengths(gtf: str):
    exons = defaultdict(list)
    gene_names = {}
    chromosomes = {}
    strands = {}

    with open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attr = fields
            attrs = parse_attrs(attr)
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_name = attrs.get("gene_name", gene_id)
            exons[gene_id].append((int(start), int(end)))
            gene_names[gene_id] = gene_name
            chromosomes[gene_id] = chrom
            strands[gene_id] = strand

    rows = []
    for gene_id in sorted(exons):
        length = merge_intervals(exons[gene_id])
        if length > 0:
            rows.append((gene_id, gene_names.get(gene_id, gene_id), chromosomes.get(gene_id, ""), strands.get(gene_id, ""), length))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build exon-union gene length table from GTF.")
    parser.add_argument("--gtf", required=True, help="Input GTF file, optionally gzipped.")
    parser.add_argument("--output", required=True, help="Output TSV file.")
    args = parser.parse_args()

    rows = build_gene_lengths(args.gtf)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write("gene_id\tgene_name\tchromosome\tstrand\tlength\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")

    if not rows:
        raise SystemExit("No exon-derived gene lengths were generated. Check the GTF file.")


if __name__ == "__main__":
    main()
