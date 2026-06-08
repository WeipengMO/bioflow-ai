#!/usr/bin/env python3
"""Collect a compact HiC-Pro QC table."""
from __future__ import annotations
import argparse
from pathlib import Path
import gzip
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hicpro-dir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--validpairs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def count_lines(path: Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt", errors="ignore") as handle:
        for _ in handle:
            n += 1
    return n


def parse_stat_files(hicpro_dir: Path, sample: str) -> dict:
    metrics = {}
    candidates = list(hicpro_dir.glob(f"**/{sample}/**/*stat*")) + list(hicpro_dir.glob(f"**/{sample}/**/*.txt"))
    key_map = {
        "valid_interaction": "valid_interaction_pairs",
        "valid pairs": "valid_interaction_pairs",
        "trans_interaction": "trans_interaction_pairs",
        "cis_interaction": "cis_interaction_pairs",
        "duplicated": "duplicated_pairs",
        "singleton": "singleton_pairs",
        "multiple": "multi_mapped_pairs",
    }
    for path in candidates:
        if not path.is_file():
            continue
        try:
            text = path.read_text(errors="ignore")
        except Exception:
            continue
        for line in text.splitlines():
            parts = re.split(r"[:\t ]+", line.strip())
            if len(parts) < 2:
                continue
            value = None
            for token in reversed(parts):
                try:
                    value = float(token.replace(",", ""))
                    break
                except ValueError:
                    continue
            if value is None:
                continue
            lower = line.lower()
            for needle, key in key_map.items():
                if needle in lower and key not in metrics:
                    metrics[key] = value
    return metrics


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    validpairs = dict(zip(samples, [Path(p) for p in args.validpairs]))
    rows = []
    for sample in samples:
        vp = validpairs.get(sample)
        metrics = parse_stat_files(Path(args.hicpro_dir), sample)
        valid_count = None
        if vp and vp.exists():
            try:
                valid_count = count_lines(vp)
            except Exception as exc:
                print(f"WARNING: failed to count valid pairs for {sample}: {exc}")
        row = {
            "sample": sample,
            "valid_pairs_file": str(vp or ""),
            "valid_pairs_counted": valid_count if valid_count is not None else "NA",
        }
        row.update(metrics)
        if "cis_interaction_pairs" in row and "trans_interaction_pairs" in row:
            trans = float(row["trans_interaction_pairs"])
            row["cis_trans_ratio"] = float(row["cis_interaction_pairs"]) / trans if trans > 0 else "Inf"
        rows.append(row)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
