#!/usr/bin/env python3
"""Identify RNaseH-sensitive regions from treatment CUT&Tag R-loop peaks."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import pyBigWig


def read_bed(path: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if start < end:
                rows.append(fields)
    return rows


def mean_signal(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> float:
    chroms = bw.chroms()
    chrom_size = chroms.get(chrom)
    if chrom_size is None:
        return 0.0
    clipped_start = max(0, min(start, chrom_size))
    clipped_end = max(0, min(end, chrom_size))
    if clipped_start >= clipped_end:
        return 0.0
    value = bw.stats(chrom, clipped_start, clipped_end, type="mean", exact=True)[0]
    if value is None or math.isnan(value):
        return 0.0
    return float(value)


def safe_ratio(treatment: float, rnaseh: float) -> float:
    if rnaseh == 0:
        return math.inf if treatment > 0 else 0.0
    return treatment / rnaseh


def format_float(value: float) -> str:
    if math.isinf(value):
        return "inf"
    return f"{value:.8g}"


def peak_name(fields: list[str], index: int) -> str:
    if len(fields) >= 4 and fields[3]:
        return fields[3]
    return f"peak_{index}"


def write_summary(path: str, sample: str, rnaseh_sample: str, peak_count: int, sensitive_count: int, args: argparse.Namespace) -> None:
    sensitive_fraction = sensitive_count / peak_count if peak_count else 0.0
    rows = [
        ("sample", sample),
        ("rnaseh_sample", rnaseh_sample),
        ("candidate_peak_count", peak_count),
        ("rnaseh_sensitive_region_count", sensitive_count),
        ("fraction_sensitive", f"{sensitive_fraction:.6f}"),
        ("scale_method", args.scale_method),
        ("min_fold_change", args.min_fold_change),
        ("min_treatment_signal", args.min_treatment_signal),
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--rnaseh-sample", required=True)
    parser.add_argument("--treatment-peaks", required=True)
    parser.add_argument("--treatment-bw", required=True)
    parser.add_argument("--rnaseh-bw", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-fold-change", type=float, default=2.0)
    parser.add_argument("--min-treatment-signal", type=float, default=0.0)
    parser.add_argument("--scale-method", required=True)
    args = parser.parse_args()

    if args.scale_method == "CPM":
        print(
            "WARNING: RNaseH-sensitive regions are being identified from CPM-normalized bigWigs. "
            "CPM is useful for browser-scale signal comparison but can obscure global RNaseH depletion. "
            "Use spikein scale when spike-in data are available."
        )

    peaks = read_bed(args.treatment_peaks)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_bed).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)

    depleted_count = 0
    fields = [
        "chrom",
        "start",
        "end",
        "name",
        "treatment_signal",
        "rnaseh_signal",
        "signal_ratio",
        "signal_difference",
        "classification",
        "scale_method",
    ]
    with pyBigWig.open(args.treatment_bw) as treatment_bw, pyBigWig.open(args.rnaseh_bw) as rnaseh_bw:
        with open(args.output_tsv, "w", newline="") as tsv_handle, open(args.output_bed, "w") as bed_handle:
            writer = csv.DictWriter(tsv_handle, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            for index, peak in enumerate(peaks, start=1):
                chrom, start_s, end_s = peak[:3]
                start = int(start_s)
                end = int(end_s)
                treatment_signal = mean_signal(treatment_bw, chrom, start, end)
                rnaseh_signal = mean_signal(rnaseh_bw, chrom, start, end)
                ratio = safe_ratio(treatment_signal, rnaseh_signal)
                diff = treatment_signal - rnaseh_signal
                sensitive = (
                    treatment_signal > rnaseh_signal
                    and treatment_signal >= args.min_treatment_signal
                    and ratio >= args.min_fold_change
                )
                classification = "rnaseh_sensitive" if sensitive else "not_sensitive"
                name = peak_name(peak, index)
                writer.writerow(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "name": name,
                        "treatment_signal": format_float(treatment_signal),
                        "rnaseh_signal": format_float(rnaseh_signal),
                        "signal_ratio": format_float(ratio),
                        "signal_difference": format_float(diff),
                        "classification": classification,
                        "scale_method": args.scale_method,
                    }
                )
                if sensitive:
                    depleted_count += 1
                    bed_handle.write(
                        "\t".join(
                            [
                                chrom,
                                str(start),
                                str(end),
                                name,
                                "0",
                                ".",
                                format_float(treatment_signal),
                                format_float(rnaseh_signal),
                                format_float(ratio),
                            ]
                        )
                        + "\n"
                    )

    write_summary(
        args.summary,
        args.sample,
        args.rnaseh_sample,
        len(peaks),
        depleted_count,
        args,
    )


if __name__ == "__main__":
    main()
