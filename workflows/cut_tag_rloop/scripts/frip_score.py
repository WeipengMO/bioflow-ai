#!/usr/bin/env python3
"""Calculate CUT&Tag FRiP with paired-end fragment semantics."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def run_count(cmd: list[str]) -> int:
    value = subprocess.check_output(cmd, text=True).strip()
    return int(value or "0")


def proper_pair_fragments(bam: str, threads: int) -> int:
    reads = run_count(["samtools", "view", "-@", str(threads), "-f", "2", "-F", "3852", "-c", bam])
    return reads // 2


def single_end_alignments(bam: str, threads: int) -> int:
    return run_count(["samtools", "view", "-@", str(threads), "-c", bam])


def write_saf_from_peaks(peaks: str, saf: str) -> int:
    count = 0
    with open(peaks) as src, open(saf, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        for line in src:
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
            if start >= end:
                continue
            count += 1
            writer.writerow([f"peak_{count}", fields[0], start + 1, end, "."])
    return count


def parse_featurecounts_assigned(path: str) -> int:
    assigned = 0
    with open(path, newline="") as handle:
        data_lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(data_lines, delimiter="\t")
    sample_columns = (reader.fieldnames or [])[6:]
    if not sample_columns:
        return 0
    sample_column = sample_columns[0]
    for row in reader:
        assigned += int(row.get(sample_column, "0") or "0")
    return assigned


def paired_end_fragments_in_peaks(bam: str, peaks: str, threads: int, workdir: str) -> int:
    saf = str(Path(workdir, "peaks.saf"))
    peak_count = write_saf_from_peaks(peaks, saf)
    if peak_count == 0:
        return 0
    output = str(Path(workdir, "featureCounts.txt"))
    subprocess.run(
        [
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
            "0",
            "-s",
            "0",
            "-p",
            "--countReadPairs",
            "-B",
            "-C",
            bam,
        ],
        check=True,
    )
    return parse_featurecounts_assigned(output)


def single_end_alignments_in_peaks(bam: str, peaks: str, threads: int) -> int:
    intersect = subprocess.Popen(["bedtools", "intersect", "-u", "-abam", bam, "-b", peaks], stdout=subprocess.PIPE)
    try:
        assert intersect.stdout is not None
        view = subprocess.run(
            ["samtools", "view", "-@", str(threads), "-c", "-"],
            stdin=intersect.stdout,
            text=True,
            stdout=subprocess.PIPE,
            check=True,
        )
        intersect.stdout.close()
        code = intersect.wait()
        if code:
            raise subprocess.CalledProcessError(code, intersect.args)
        return int(view.stdout.strip() or "0")
    finally:
        if intersect.poll() is None:
            intersect.kill()


def write_frip(path: str, sample: str, total: int, in_peaks: int, counting_mode: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frip = in_peaks / total if total else 0
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "total_fragments", "fragments_in_peaks", "frip", "counting_mode"])
        writer.writerow([sample, total, in_peaks, f"{frip:.6f}", counting_mode])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=["pe", "se"], required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    if args.mode == "pe":
        total = proper_pair_fragments(args.bam, args.threads)
        counting_mode = "paired_end_fragments"
        if Path(args.peaks).exists() and Path(args.peaks).stat().st_size:
            with tempfile.TemporaryDirectory(prefix=f"{args.sample}.frip.") as tmpdir:
                in_peaks = paired_end_fragments_in_peaks(args.bam, args.peaks, args.threads, tmpdir)
        else:
            in_peaks = 0
    else:
        total = single_end_alignments(args.bam, args.threads)
        counting_mode = "single_end_alignments"
        in_peaks = (
            single_end_alignments_in_peaks(args.bam, args.peaks, args.threads)
            if Path(args.peaks).exists() and Path(args.peaks).stat().st_size
            else 0
        )
    write_frip(args.output, args.sample, total, in_peaks, counting_mode)


if __name__ == "__main__":
    main()
