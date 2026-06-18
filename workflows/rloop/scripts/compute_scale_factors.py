#!/usr/bin/env python3
"""Compute sample-level bigWig scale factors from workflow BAMs."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def count_alignments(bam: str) -> int:
    value = subprocess.check_output(["samtools", "view", "-c", bam], text=True).strip()
    return int(value)


def load_source_tsv(path: str) -> dict[str, float]:
    if not path:
        raise ValueError("--source-tsv is required when --method tsv")
    factors: dict[str, float] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample", "scale_factor"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError("source TSV is missing column(s): " + ", ".join(sorted(missing)))
        for row in reader:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            factor = float(row["scale_factor"])
            if factor <= 0:
                raise ValueError(f"scale_factor must be > 0 for sample {sample}")
            factors[sample] = factor
    return factors


def load_spikein_counts(path: str) -> dict[str, int]:
    if not path:
        raise ValueError("--spikein-counts is required when --method spikein")
    counts: dict[str, int] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample", "spikein_aligned"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError("spikein counts TSV is missing column(s): " + ", ".join(sorted(missing)))
        for row in reader:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            count = int(row["spikein_aligned"])
            if count < 0:
                raise ValueError(f"spikein_aligned must be >= 0 for sample {sample}")
            counts[sample] = count
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Comma-separated sample names in BAM order.")
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--mode", choices=["pe", "se"], required=True)
    parser.add_argument("--method", choices=["effective_fragments", "raw", "tsv", "spikein"], required=True)
    parser.add_argument("--source-tsv", default="")
    parser.add_argument("--spikein-counts", default="", help="TSV with sample and spikein_aligned columns (for spikein method).")
    parser.add_argument("--spikein-method", choices=["ratio", "rpm"], default="ratio")
    parser.add_argument("--spikein-reference-sample", default="", help="Sample whose spikein count defines scale=1.0 (ratio mode).")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    samples = [sample for sample in args.samples.split(",") if sample]
    if len(samples) != len(args.bams):
        raise ValueError(f"Expected {len(samples)} BAMs for samples but got {len(args.bams)}")

    source_factors = load_source_tsv(args.source_tsv) if args.method == "tsv" else {}
    spikein_counts: dict[str, int] = {}
    spikein_ref_count = 0
    if args.method == "spikein":
        spikein_counts = load_spikein_counts(args.spikein_counts)
        if args.spikein_method == "ratio":
            if not args.spikein_reference_sample:
                raise ValueError("--spikein-reference-sample is required for ratio mode")
            if args.spikein_reference_sample not in spikein_counts:
                raise ValueError(f"spikein reference sample not found in counts: {args.spikein_reference_sample}")
            spikein_ref_count = spikein_counts[args.spikein_reference_sample]
            if spikein_ref_count <= 0:
                raise ValueError(f"spikein reference sample has 0 alignments: {args.spikein_reference_sample}")
    rows = []
    for sample, bam in zip(samples, args.bams):
        alignments = count_alignments(bam)
        effective_fragments = alignments / 2 if args.mode == "pe" else alignments
        if args.method == "raw":
            scale_factor = 1.0
        elif args.method == "tsv":
            if sample not in source_factors:
                raise ValueError(f"source TSV is missing sample: {sample}")
            scale_factor = source_factors[sample]
        elif args.method == "spikein":
            if sample not in spikein_counts:
                raise ValueError(f"spikein counts is missing sample: {sample}")
            count = spikein_counts[sample]
            if count <= 0:
                raise ValueError(f"spikein alignment count is 0 for sample: {sample}")
            if args.spikein_method == "ratio":
                scale_factor = spikein_ref_count / count
            else:
                scale_factor = 1_000_000.0 / count
        else:
            if effective_fragments <= 0:
                raise ValueError(f"No effective fragments for sample: {sample}")
            scale_factor = 1_000_000.0 / effective_fragments
        rows.append(
            {
                "sample": sample,
                "alignments": int(alignments),
                "scale_factor": f"{scale_factor:.12g}",
                "effective_fragments": f"{effective_fragments:.6f}",
                "method": args.method,
                "bam": bam,
            }
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        fields = ["sample", "alignments", "scale_factor", "effective_fragments", "method", "bam"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
