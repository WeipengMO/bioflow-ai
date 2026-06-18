#!/usr/bin/env python3
"""Build a shared peak universe and quantify it with featureCounts."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from pathlib import Path


def run(cmd: list[str], stdout=None) -> None:
    subprocess.run(cmd, check=True, stdout=stdout)


def read_bed3(path: str):
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
                yield fields[0], start, end


def build_universe(
    peak_files: list[str],
    output: str,
    min_width: int,
    merge_distance: int,
    universe_mode: str,
    min_support: int,
) -> None:
    intervals = []
    for sample_index, path in enumerate(peak_files):
        for chrom, start, end in read_bed3(path):
            intervals.append((chrom, start, end, 1 << sample_index))
    intervals.sort(key=lambda row: (row[0], row[1], row[2]))
    merged: list[list[object]] = []
    for chrom, start, end, sample_bit in intervals:
        if not merged or merged[-1][0] != chrom or start > int(merged[-1][2]) + merge_distance:
            merged.append([chrom, start, end, sample_bit])
        else:
            merged[-1][2] = max(int(merged[-1][2]), end)
            merged[-1][3] = int(merged[-1][3]) | sample_bit

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        idx = 0
        for chrom, start, end, sample_bits in merged:
            if int(end) - int(start) < min_width:
                continue
            support = int(sample_bits).bit_count()
            if universe_mode == "consensus":
                if support < min_support:
                    continue
            idx += 1
            handle.write(f"{chrom}\t{start}\t{end}\tpeak_{idx}\t{support}\t.\n")


def write_annotation_input(universe: str, output: str) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(universe) as src, open(output, "w") as out:
        for idx, line in enumerate(src, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                name = fields[3] if len(fields) >= 4 and fields[3] else f"peak_{idx}"
                out.write("\t".join([fields[0], fields[1], fields[2], name]) + "\n")


def write_saf(universe: str, output: str) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(universe) as src, open(output, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for idx, line in enumerate(src, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name = fields[3] if len(fields) >= 4 and fields[3] else f"peak_{idx}"
            start_1based = int(fields[1]) + 1
            out.write("\t".join([name, fields[0], str(start_1based), fields[2], "."]) + "\n")


def read_universe(universe: str) -> list[list[str]]:
    rows = []
    with open(universe) as handle:
        for idx, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                rows.append([fields[0], fields[1], fields[2], fields[3] if len(fields) >= 4 else f"peak_{idx}"])
    return rows


def run_featurecounts(
    saf: str,
    bams: list[str],
    output: str,
    assay_mode: str,
    min_mapq: int,
    threads: int,
    extra: str,
) -> None:
    cmd = [
        "featureCounts",
        "-T",
        str(max(1, threads)),
        "-F",
        "SAF",
        "-a",
        saf,
        "-o",
        output,
        "-Q",
        str(min_mapq),
        "-s",
        "0",
    ]
    if assay_mode == "pe":
        cmd.extend(["-p", "--countReadPairs", "-B", "-C"])
    if extra:
        cmd.extend(shlex.split(extra))
    cmd.extend(bams)
    subprocess.run(cmd, check=True)


def parse_featurecounts(path: str, samples: list[str], universe_rows: list[list[str]]) -> list[list[int]]:
    with open(path, newline="") as handle:
        data_lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(data_lines, delimiter="\t")
    sample_columns = (reader.fieldnames or [])[6:]
    if len(sample_columns) != len(samples):
        raise ValueError(
            f"featureCounts output has {len(sample_columns)} sample columns, expected {len(samples)}."
        )
    rows_by_peak = {}
    for row in reader:
        rows_by_peak[row["Geneid"]] = [int(row[column]) for column in sample_columns]
    matrix = []
    for row in universe_rows:
        matrix.append(rows_by_peak.get(row[3], [0] * len(samples)))
    return matrix


def read_spikein_scale_factors(path: str | None) -> dict[str, float]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    values = {}
    with p.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row.get("sample", "")
            scale = row.get("spikein_scale_factor", "")
            if sample and scale and scale != "NA":
                values[sample] = float(scale)
    return values


def write_matrix(path: str, universe_rows: list[list[str]], samples: list[str], matrix: list[list[float]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["peak_id", "chrom", "start", "end", *samples])
        for row, values in zip(universe_rows, matrix):
            writer.writerow([row[3], row[0], row[1], row[2], *values])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--mode", choices=["pe", "se"], required=True)
    parser.add_argument("--peak-files", nargs="+", required=True)
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--universe-mode", choices=["union", "consensus"], default="consensus")
    parser.add_argument("--min-support", type=int, default=2)
    parser.add_argument("--min-peak-width", type=int, default=50)
    parser.add_argument("--merge-distance", type=int, default=100)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--featurecounts-extra", nargs="?", const="", default="")
    parser.add_argument("--normalization-metrics", default="")
    parser.add_argument("--output-universe", required=True)
    parser.add_argument("--output-saf", required=True)
    parser.add_argument("--output-featurecounts", required=True)
    parser.add_argument("--output-raw", required=True)
    parser.add_argument("--output-cpm", required=True)
    parser.add_argument("--output-spikein", required=True)
    parser.add_argument("--output-annotation-bed", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    samples = [sample for sample in args.samples.split(",") if sample]
    if len(samples) != len(args.bams):
        raise ValueError(f"Expected {len(samples)} BAMs, got {len(args.bams)}")

    build_universe(
        args.peak_files,
        args.output_universe,
        args.min_peak_width,
        args.merge_distance,
        args.universe_mode,
        args.min_support,
    )
    write_annotation_input(args.output_universe, args.output_annotation_bed)
    write_saf(args.output_universe, args.output_saf)
    universe_rows = read_universe(args.output_universe)

    run_featurecounts(
        args.output_saf,
        args.bams,
        args.output_featurecounts,
        args.mode,
        args.min_mapq,
        args.threads,
        args.featurecounts_extra,
    )
    raw_matrix = parse_featurecounts(args.output_featurecounts, samples, universe_rows)
    cpm_matrix = []
    totals = [sum(row[index] for row in raw_matrix) for index in range(len(samples))]
    for row in raw_matrix:
        cpm_matrix.append([f"{(value / total * 1_000_000) if total else 0:.8g}" for value, total in zip(row, totals)])

    write_matrix(args.output_raw, universe_rows, samples, raw_matrix)
    write_matrix(args.output_cpm, universe_rows, samples, cpm_matrix)

    scale_factors = read_spikein_scale_factors(args.normalization_metrics)
    spikein_matrix = []
    for row in raw_matrix:
        spikein_matrix.append([f"{value * scale_factors.get(sample, 0.0):.8g}" for value, sample in zip(row, samples)])
    write_matrix(args.output_spikein, universe_rows, samples, spikein_matrix)


if __name__ == "__main__":
    main()
