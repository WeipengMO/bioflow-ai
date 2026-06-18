#!/usr/bin/env python3
"""Collect generic spike-in normalization metrics for CUT&Tag R-loop samples."""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
from pathlib import Path


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def count_fastq_reads(paths: list[str]) -> int:
    total = 0
    for path in paths:
        if not path:
            continue
        lines = 0
        with open_maybe_gzip(path) as handle:
            for lines, _ in enumerate(handle, start=1):
                pass
        total += lines // 4
    return total


def samtools_count(args: list[str]) -> int:
    value = subprocess.check_output(["samtools", "view", "-c", *args], text=True).strip()
    return int(value or "0")


def mapped_reads(bam: str) -> int:
    return samtools_count(["-F", "4", bam])


def proper_pair_fragments(bam: str) -> int:
    return samtools_count(["-f", "2", "-F", "3852", bam]) // 2


def unique_fragments(bam: str, mode: str) -> int:
    if mode == "pe":
        return samtools_count(["-f", "2", "-F", "3852", bam]) // 2
    return samtools_count(["-F", "3844", bam])


def read_fastq_manifest(path: str) -> dict[str, list[str]]:
    fastqs: dict[str, list[str]] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample", "read1"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError("FASTQ metrics manifest is missing column(s): " + ", ".join(sorted(missing)))
        for row in reader:
            sample = row["sample"].strip()
            if not sample:
                continue
            reads = [row["read1"].strip()]
            read2 = (row.get("read2") or "").strip()
            if read2:
                reads.append(read2)
            fastqs[sample] = reads
    return fastqs


def read_spikein_groups(path: str) -> dict[str, dict[str, str]]:
    groups: dict[str, dict[str, str]] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "spikein_group", "spikein_reference_sample"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Spike-in group manifest is missing column(s): " + ", ".join(sorted(missing)))
        for row in reader:
            sample = row["sample"].strip()
            if not sample:
                continue
            groups[sample] = {
                "spikein_group": row["spikein_group"].strip(),
                "spikein_reference_sample": row["spikein_reference_sample"].strip(),
            }
    return groups


def format_float(value: float) -> str:
    return f"{value:.8g}"


def write_warnings(path: str, text_path: str, rows: list[dict[str, str]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "spikein_genome",
        "spikein_count",
        "spikein_fraction",
        "threshold",
        "severity",
        "message",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with open(text_path, "w") as handle:
        if rows:
            for row in rows:
                handle.write(
                    "{severity}: {sample}: {message} "
                    "(spikein_genome={spikein_genome}, spikein_count={spikein_count}, "
                    "spikein_fraction={spikein_fraction}, threshold={threshold})\n".format(
                        **row
                    )
                )
        else:
            handle.write("No spike-in warnings.\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Comma-separated sample names in BAM order.")
    parser.add_argument("--mode", choices=["pe", "se"], required=True)
    parser.add_argument("--fastq-manifest", required=True)
    parser.add_argument("--spikein-groups", required=True)
    parser.add_argument("--human-bams", nargs="+", required=True)
    parser.add_argument("--human-unique-bams", nargs="+", required=True)
    parser.add_argument("--spikein-bams", nargs="+", required=True)
    parser.add_argument("--spikein-genome", default="ecoli")
    parser.add_argument("--spikein-counting-mode", choices=["mapped_reads", "proper_pair_fragments"], default="mapped_reads")
    parser.add_argument("--min-spikein-reads", type=int, default=1000)
    parser.add_argument("--warn-low-fraction", type=float, default=0.001)
    parser.add_argument("--output", required=True)
    parser.add_argument("--warning-tsv", required=True)
    parser.add_argument("--warning-txt", required=True)
    args = parser.parse_args()

    samples = [sample for sample in args.samples.split(",") if sample]
    for label, values in {
        "human-bams": args.human_bams,
        "human-unique-bams": args.human_unique_bams,
        "spikein-bams": args.spikein_bams,
    }.items():
        if len(values) != len(samples):
            raise ValueError(f"Expected {len(samples)} {label} values, got {len(values)}")
    fastqs = read_fastq_manifest(args.fastq_manifest)
    missing_fastqs = sorted(set(samples) - set(fastqs))
    if missing_fastqs:
        raise ValueError("FASTQ metrics manifest is missing sample(s): " + ", ".join(missing_fastqs))
    spikein_groups = read_spikein_groups(args.spikein_groups)
    missing_groups = sorted(set(samples) - set(spikein_groups))
    if missing_groups:
        raise ValueError("Spike-in group manifest is missing sample(s): " + ", ".join(missing_groups))

    raw_rows = []
    warnings = []
    for sample, human_bam, human_unique_bam, spikein_bam in zip(
        samples, args.human_bams, args.human_unique_bams, args.spikein_bams
    ):
        total_reads = count_fastq_reads(fastqs[sample])
        human_mapped = mapped_reads(human_bam)
        human_unique = unique_fragments(human_unique_bam, args.mode)
        spikein_mapped = mapped_reads(spikein_bam)
        spikein_proper_pairs = proper_pair_fragments(spikein_bam) if args.mode == "pe" else 0
        spikein_count = spikein_proper_pairs if args.spikein_counting_mode == "proper_pair_fragments" and args.mode == "pe" else spikein_mapped
        spikein_fraction = spikein_mapped / total_reads if total_reads else 0.0
        group_info = spikein_groups[sample]

        warning_level = "ok"
        warning_message = ""
        if spikein_count == 0:
            warning_level = "error"
            warning_message = "Spike-in count is 0; spike-in scaling is invalid."
        elif spikein_count < args.min_spikein_reads or spikein_fraction < args.warn_low_fraction:
            warning_level = "warning"
            warning_message = "Spike-in support is below the configured count or fraction threshold."

        if warning_level != "ok":
            warnings.append(
                {
                    "sample": sample,
                    "spikein_genome": args.spikein_genome,
                    "spikein_count": str(spikein_count),
                    "spikein_fraction": format_float(spikein_fraction),
                    "threshold": f"min_spikein_reads={args.min_spikein_reads};warn_low_fraction={args.warn_low_fraction}",
                    "severity": warning_level,
                    "message": warning_message,
                }
            )

        raw_rows.append(
            {
                "sample": sample,
                "human_total_reads": total_reads,
                "human_mapped_reads": human_mapped,
                "human_unique_fragments": human_unique,
                "spikein_genome": args.spikein_genome,
                "spikein_mapped_reads": spikein_mapped,
                "spikein_proper_pair_fragments": spikein_proper_pairs if args.mode == "pe" else "NA",
                "spikein_count": spikein_count,
                "spikein_fraction": format_float(spikein_fraction),
                "spikein_counting_mode": args.spikein_counting_mode,
                "spikein_group": group_info["spikein_group"],
                "spikein_reference_sample": group_info["spikein_reference_sample"],
                "spikein_reference_count": "NA",
                "spikein_reference_human_unique_fragments": "NA",
                "spikein_raw_scale_factor": "NA",
                "spikein_unit_scale_factor": "NA",
                "spikein_scale_factor": "NA",
                "spikein_warning": warning_level if warning_level != "ok" else "",
                "warning_message": warning_message,
            }
        )

    rows_by_sample = {str(row["sample"]): row for row in raw_rows}
    reported_reference_errors = set()
    for row in raw_rows:
        reference_sample = str(row["spikein_reference_sample"])
        reference = rows_by_sample.get(reference_sample)
        if reference is None:
            raise ValueError(f"Spike-in reference sample is not declared: {reference_sample}")
        reference_count = int(reference["spikein_count"])
        reference_unique = int(reference["human_unique_fragments"])
        row["spikein_reference_count"] = reference_count
        row["spikein_reference_human_unique_fragments"] = reference_unique
        if reference_count <= 0 and (reference_sample, "spikein") not in reported_reference_errors:
            reported_reference_errors.add((reference_sample, "spikein"))
            warnings.append(
                {
                    "sample": reference_sample,
                    "spikein_genome": args.spikein_genome,
                    "spikein_count": str(reference_count),
                    "spikein_fraction": str(reference["spikein_fraction"]),
                    "threshold": "reference_spikein_count>0",
                    "severity": "error",
                    "message": "Spike-in reference sample has no usable spike-in reads/fragments.",
                }
            )
        if reference_unique <= 0 and (reference_sample, "human_unique") not in reported_reference_errors:
            reported_reference_errors.add((reference_sample, "human_unique"))
            warnings.append(
                {
                    "sample": reference_sample,
                    "spikein_genome": args.spikein_genome,
                    "spikein_count": str(reference_count),
                    "spikein_fraction": str(reference["spikein_fraction"]),
                    "threshold": "reference_human_unique_fragments>0",
                    "severity": "error",
                    "message": "Spike-in reference sample has no usable human unique fragments for unit scaling.",
                }
            )
        count = int(row["spikein_count"])
        if count > 0 and reference_count > 0 and reference_unique > 0:
            raw_scale = reference_count / count
            unit_scale = 1_000_000 / reference_unique
            row["spikein_raw_scale_factor"] = format_float(raw_scale)
            row["spikein_unit_scale_factor"] = format_float(unit_scale)
            row["spikein_scale_factor"] = format_float(raw_scale * unit_scale)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "human_total_reads",
        "human_mapped_reads",
        "human_unique_fragments",
        "spikein_genome",
        "spikein_mapped_reads",
        "spikein_proper_pair_fragments",
        "spikein_count",
        "spikein_fraction",
        "spikein_counting_mode",
        "spikein_group",
        "spikein_reference_sample",
        "spikein_reference_count",
        "spikein_reference_human_unique_fragments",
        "spikein_raw_scale_factor",
        "spikein_unit_scale_factor",
        "spikein_scale_factor",
        "spikein_warning",
        "warning_message",
    ]
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(raw_rows)

    write_warnings(args.warning_tsv, args.warning_txt, warnings)
    if any(row["severity"] == "error" for row in warnings):
        raise SystemExit("Spike-in normalization has error-level warnings; see warning files.")


if __name__ == "__main__":
    main()
