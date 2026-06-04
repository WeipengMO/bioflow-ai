#!/usr/bin/env python3
"""Filter gffcompare-annotated transcripts into candidate lncRNAs."""

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
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr)}


def format_attrs(attrs: Dict[str, str]) -> str:
    return " ".join(f'{key} "{value}";' for key, value in attrs.items())


def read_fasta_seqnames(path: str | None) -> set[str] | None:
    if not path:
        return None
    seqnames = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                seqnames.add(line[1:].split()[0])
    if not seqnames:
        raise ValueError(f"No FASTA records found in {path}")
    return seqnames


def transcript_id(fields: List[str]) -> str | None:
    attrs = parse_attrs(fields[8])
    return attrs.get("transcript_id")


def exon_stats(records: Iterable[List[str]]) -> Tuple[int, int]:
    exons = []
    for fields in records:
        if fields[2] == "exon":
            exons.append((int(fields[3]), int(fields[4])))
    if not exons:
        return 0, 0
    exons = sorted(exons)
    merged = []
    cur_start, cur_end = exons[0]
    for start, end in exons[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    length = sum(end - start + 1 for start, end in merged)
    return length, len(exons)


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter candidate lncRNAs from gffcompare annotated GTF.")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output-gtf", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-length", type=int, default=200)
    parser.add_argument("--min-exons", type=int, default=2)
    parser.add_argument("--class-codes", default="u,i,x,o,e")
    parser.add_argument("--fasta", default=None, help="Optional genome FASTA. Transcripts on missing seqnames are skipped.")
    args = parser.parse_args()

    keep_codes = {code.strip() for code in args.class_codes.split(",") if code.strip()}
    allowed_seqnames = read_fasta_seqnames(args.fasta)
    records_by_tx: Dict[str, List[List[str]]] = defaultdict(list)
    attrs_by_tx: Dict[str, Dict[str, str]] = {}

    with open_text(args.gtf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            tx_id = transcript_id(fields)
            if not tx_id:
                continue
            records_by_tx[tx_id].append(fields)
            attrs_by_tx.setdefault(tx_id, parse_attrs(fields[8]))

    kept = []
    for tx_id, records in sorted(records_by_tx.items()):
        attrs = attrs_by_tx.get(tx_id, {})
        class_code = attrs.get("class_code", "")
        length, exon_count = exon_stats(records)
        if allowed_seqnames is not None and any(fields[0] not in allowed_seqnames for fields in records):
            continue
        if class_code in keep_codes and length >= args.min_length and exon_count >= args.min_exons:
            kept.append((tx_id, class_code, length, exon_count, records))

    Path(args.output_gtf).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_gtf, "w") as out_gtf:
        for tx_id, class_code, length, exon_count, records in kept:
            for fields in records:
                attrs = parse_attrs(fields[8])
                attrs.setdefault("gene_name", attrs.get("gene_id", tx_id))
                attrs["lncRNA_candidate"] = "yes"
                attrs["class_code"] = class_code
                attrs["transcript_length"] = str(length)
                attrs["exon_count"] = str(exon_count)
                fields = fields.copy()
                fields[8] = format_attrs(attrs)
                out_gtf.write("\t".join(fields) + "\n")

    with open(args.summary, "w") as out_summary:
        out_summary.write("transcript_id\tgene_id\tclass_code\tlength\texon_count\n")
        for tx_id, class_code, length, exon_count, records in kept:
            attrs = attrs_by_tx.get(tx_id, {})
            out_summary.write(f"{tx_id}\t{attrs.get('gene_id', '')}\t{class_code}\t{length}\t{exon_count}\n")

    if not kept:
        raise SystemExit("No candidate lncRNAs passed filters. Relax lncrna_filter settings or inspect gffcompare output.")


if __name__ == "__main__":
    main()
