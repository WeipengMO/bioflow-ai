#!/usr/bin/env python3
"""Build a consensus peak universe with sample support metadata."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(order=True)
class Interval:
    chrom: str
    start: int
    end: int
    sample: str


@dataclass
class MergedPeak:
    chrom: str
    start: int
    end: int
    samples: set[str] = field(default_factory=set)

    def add(self, interval: Interval) -> None:
        self.end = max(self.end, interval.end)
        self.samples.add(interval.sample)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge BED-like peak files into a consensus peak universe."
    )
    parser.add_argument(
        "--sample",
        nargs=2,
        action="append",
        metavar=("SAMPLE", "PEAK_FILE"),
        required=True,
        help="Sample name and BED/narrowPeak/broadPeak file. Repeat once per sample.",
    )
    parser.add_argument("--bed", required=True, help="Output BED6 path.")
    parser.add_argument("--tsv", required=True, help="Output TSV path with support metadata.")
    parser.add_argument("--matrix", required=True, help="Output sample-by-peak support matrix path.")
    parser.add_argument("--saf", required=True, help="Output SAF path for read counting.")
    parser.add_argument("--min-support-count", type=int, default=1, help="Minimum supporting sample count.")
    parser.add_argument(
        "--min-support-fraction",
        type=float,
        default=0.0,
        help="Minimum supporting sample fraction in [0, 1].",
    )
    return parser.parse_args()


def validate_thresholds(min_support_count: int, min_support_fraction: float) -> None:
    if min_support_count < 1:
        raise ValueError("--min-support-count must be at least 1.")
    if min_support_fraction < 0 or min_support_fraction > 1:
        raise ValueError("--min-support-fraction must be between 0 and 1.")


def read_intervals(sample: str, peak_file: str) -> list[Interval]:
    intervals = []
    with Path(peak_file).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{peak_file}:{line_number} has non-integer coordinates.") from exc
            if start < 0 or start >= end:
                continue
            intervals.append(Interval(fields[0], start, end, sample))
    return intervals


def merge_intervals(intervals: list[Interval]) -> list[MergedPeak]:
    merged: list[MergedPeak] = []
    for interval in sorted(intervals):
        if not merged or interval.chrom != merged[-1].chrom or interval.start > merged[-1].end:
            merged.append(MergedPeak(interval.chrom, interval.start, interval.end, {interval.sample}))
        else:
            merged[-1].add(interval)
    return merged


def filter_peaks(
    peaks: list[MergedPeak],
    sample_count: int,
    min_support_count: int,
    min_support_fraction: float,
) -> list[MergedPeak]:
    required_count = max(min_support_count, int(min_support_fraction * sample_count + 0.999999))
    return [peak for peak in peaks if len(peak.samples) >= required_count]


def write_outputs(
    peaks: list[MergedPeak],
    samples: list[str],
    bed_path: str,
    tsv_path: str,
    matrix_path: str,
    saf_path: str,
) -> None:
    for output_path in [bed_path, tsv_path, matrix_path, saf_path]:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with (
        Path(bed_path).open("w", newline="") as bed,
        Path(tsv_path).open("w", newline="") as tsv,
        Path(matrix_path).open("w", newline="") as matrix,
        Path(saf_path).open("w", newline="") as saf,
    ):
        tsv_writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        matrix_writer = csv.writer(matrix, delimiter="\t", lineterminator="\n")
        saf_writer = csv.writer(saf, delimiter="\t", lineterminator="\n")

        tsv_writer.writerow(["peak_id", "chrom", "start", "end", "support_count", "support_fraction", "support_samples"])
        matrix_writer.writerow(["peak_id", "chrom", "start", "end", *samples])
        saf_writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])

        for index, peak in enumerate(peaks, start=1):
            peak_id = f"consensus_peak_{index}"
            support_samples = sorted(peak.samples)
            support_count = len(support_samples)
            support_fraction = support_count / len(samples)
            bed.write(
                "\t".join([peak.chrom, str(peak.start), str(peak.end), peak_id, str(support_count), "."])
                + "\n"
            )
            tsv_writer.writerow(
                [
                    peak_id,
                    peak.chrom,
                    peak.start,
                    peak.end,
                    support_count,
                    f"{support_fraction:.6g}",
                    ",".join(support_samples),
                ]
            )
            matrix_writer.writerow(
                [
                    peak_id,
                    peak.chrom,
                    peak.start,
                    peak.end,
                    *["1" if sample in peak.samples else "0" for sample in samples],
                ]
            )
            saf_writer.writerow([peak_id, peak.chrom, peak.start + 1, peak.end, "."])


def main() -> None:
    args = parse_args()
    validate_thresholds(args.min_support_count, args.min_support_fraction)
    samples = [sample for sample, _ in args.sample]
    if len(samples) != len(set(samples)):
        raise ValueError("Sample names must be unique.")

    intervals = []
    for sample, peak_file in args.sample:
        intervals.extend(read_intervals(sample, peak_file))
    peaks = filter_peaks(
        merge_intervals(intervals),
        len(samples),
        args.min_support_count,
        args.min_support_fraction,
    )
    if not peaks:
        raise ValueError("No consensus peaks passed the support thresholds.")
    write_outputs(peaks, samples, args.bed, args.tsv, args.matrix, args.saf)
    print(f"Wrote {len(peaks)} consensus peaks from {len(samples)} sample(s).")


if __name__ == "__main__":
    main()
