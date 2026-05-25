#!/usr/bin/env python3
"""Create a FASTQ manifest from a raw_data directory using common read-pair suffixes."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIRED_FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:[._-]R?|[._-])(?P<read>[12])(?:_001)?(?P<suffix>\.f(?:ast)?q(?:\.gz)?)$",
    re.IGNORECASE,
)


def discover(raw_dir: Path) -> dict[str, dict[str, str]]:
    manifest: dict[str, dict[str, str]] = {}
    for path in sorted(p for p in raw_dir.iterdir() if p.is_file()):
        if not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        match = PAIRED_FASTQ_RE.match(path.name)
        if not match:
            continue
        sample = match.group("sample")
        read = "read1" if match.group("read") == "1" else "read2"
        manifest.setdefault(sample, {})[read] = str(path)
    return manifest


def write_manifest(path: Path, manifest: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "read1", "read2"], delimiter="\t")
        writer.writeheader()
        for sample, reads in sorted(manifest.items()):
            writer.writerow({"sample": sample, "read1": reads.get("read1", ""), "read2": reads.get("read2", "")})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-data-dir", default="raw_data")
    parser.add_argument("--output", default="config/fastq_manifest.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_dir = Path(args.raw_data_dir)
    if not raw_dir.exists():
        raise SystemExit(f"raw data directory does not exist: {raw_dir}")
    manifest = discover(raw_dir)
    write_manifest(Path(args.output), manifest)
    print(f"Wrote {len(manifest)} samples to {args.output}")


if __name__ == "__main__":
    main()
