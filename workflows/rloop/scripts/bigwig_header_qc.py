#!/usr/bin/env python3
"""Write compact bigWig header metadata for R-loop signal-track QC."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pyBigWig


def bigwig_header(path: str) -> dict[str, object]:
    with pyBigWig.open(path) as bw:
        header = bw.header()
    return header or {}


def write_row(writer: csv.DictWriter, sample: str, track_type: str, path: str, scale_factor: str, normalization: str) -> None:
    header = bigwig_header(path)
    writer.writerow(
        {
            "sample": sample,
            "track_type": track_type,
            "path": path,
            "version": header.get("version", ""),
            "nLevels": header.get("nLevels", ""),
            "nBasesCovered": header.get("nBasesCovered", ""),
            "minVal": header.get("minVal", ""),
            "maxVal": header.get("maxVal", ""),
            "sumData": header.get("sumData", ""),
            "sumSquared": header.get("sumSquared", ""),
            "scale_factor": scale_factor,
            "normalization": normalization,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--cpm", required=True)
    parser.add_argument("--common-scale", default="")
    parser.add_argument("--scale-factor", default="")
    parser.add_argument("--common-scale-normalization", default="")
    parser.add_argument("--debug-raw", default="")
    parser.add_argument("--debug-dedup", default="")
    parser.add_argument("--dedup-scale-factor", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "track_type",
        "path",
        "version",
        "nLevels",
        "nBasesCovered",
        "minVal",
        "maxVal",
        "sumData",
        "sumSquared",
        "scale_factor",
        "normalization",
    ]
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        write_row(writer, args.sample, "CPM", args.cpm, "", "CPM")
        if args.common_scale:
            write_row(
                writer,
                args.sample,
                "common_scale",
                args.common_scale,
                args.scale_factor,
                args.common_scale_normalization,
            )
        if args.debug_raw:
            write_row(writer, args.sample, "debug_raw_scale", args.debug_raw, "1", "None")
        if args.debug_dedup:
            write_row(
                writer,
                args.sample,
                "debug_dedup_effective_fragments",
                args.debug_dedup,
                args.dedup_scale_factor,
                "None",
            )


if __name__ == "__main__":
    main()
