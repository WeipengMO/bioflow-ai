#!/usr/bin/env python3
"""Generate a generic ChIPQC sample sheet from BioFlowAI ChIP-seq metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


COLUMNS = [
    "SampleID",
    "Tissue",
    "Factor",
    "Condition",
    "Treatment",
    "Replicate",
    "bamReads",
    "ControlID",
    "bamControl",
    "Peaks",
    "PeakCaller",
]
MISSING_VALUE = "NA"


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        return {}
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(item) for item in value]


def unique_list(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def control_map(samples_config: dict[str, Any]) -> tuple[dict[str, str], list[str]]:
    treatments = unique_list(as_list(samples_config.get("treatments", [])))
    controls = samples_config.get("controls", {}) or {}
    strategy = controls.get("strategy", "none")

    if strategy == "none":
        return {}, treatments

    if strategy == "pooled":
        pooled = controls.get("pooled", {}) or {}
        if not pooled:
            raise ValueError("controls.strategy is 'pooled', but controls.pooled is empty")

        mapping: dict[str, str] = {}
        treatment_set = set(treatments)
        for pool_name, pool_config in pooled.items():
            pool_treatments = as_list((pool_config or {}).get("treatments", []))
            pool_controls = as_list((pool_config or {}).get("controls", []))
            if not pool_treatments:
                raise ValueError(f"controls.pooled.{pool_name}.treatments is empty")
            if not pool_controls:
                raise ValueError(f"controls.pooled.{pool_name}.controls is empty")
            unknown = sorted(set(pool_treatments) - treatment_set)
            if unknown:
                raise ValueError(
                    f"controls.pooled.{pool_name}.treatments contains unknown treatment(s): "
                    + ", ".join(unknown)
                )
            duplicates = sorted(set(pool_treatments) & set(mapping))
            if duplicates:
                raise ValueError(
                    "Treatment(s) assigned to more than one pooled control: "
                    + ", ".join(duplicates)
                )
            for treatment in pool_treatments:
                mapping[treatment] = str(pool_name)

        missing = sorted(treatment_set - set(mapping))
        if missing:
            raise ValueError(
                "controls.strategy is 'pooled', but no pooled control was assigned to treatment(s): "
                + ", ".join(missing)
            )
        return mapping, treatments

    if strategy == "matched":
        matched = {str(key): str(value) for key, value in (controls.get("matched", {}) or {}).items()}
        missing = sorted(set(treatments) - set(matched))
        if missing:
            raise ValueError("controls.matched is missing treatment(s): " + ", ".join(missing))
        unknown = sorted(set(matched) - set(treatments))
        if unknown:
            raise ValueError("controls.matched contains unknown treatment(s): " + ", ".join(unknown))
        return matched, treatments

    raise ValueError("controls.strategy must be one of: none, pooled, matched")


def replicate_map(replicate_config: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for samples in replicate_config.values():
        for index, sample in enumerate(as_list(samples), start=1):
            mapping[str(sample)] = str(index)
    return mapping


def peak_path(sample: str, peak_mode: str, peak_type: str) -> str:
    if peak_type not in {"narrow", "broad"}:
        raise ValueError("--peak-type must be 'narrow' or 'broad'")
    if peak_mode not in {"with_control", "without_control"}:
        raise ValueError("--peak-mode must be 'with_control' or 'without_control'")

    suffix = "_no_control" if peak_mode == "without_control" else ""
    extension = "narrowPeak" if peak_type == "narrow" else "broadPeak"
    return f"macs3_results/{peak_type}{suffix}/{sample}_peaks.{extension}"


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    samples_config = load_yaml(args.samples)
    replicate_config = load_yaml(args.replicates)
    controls, treatments = control_map(samples_config)
    replicates = replicate_map(replicate_config)

    rows: list[dict[str, str]] = []
    for treatment in treatments:
        control_id = controls.get(treatment, args.no_control_value)
        bam_control = (
            f"{args.bam_dir}/{control_id}.sorted.rmdup.bam"
            if control_id != args.no_control_value
            else args.no_control_value
        )

        rows.append(
            {
                "SampleID": treatment,
                "Tissue": MISSING_VALUE,
                "Factor": MISSING_VALUE,
                "Condition": MISSING_VALUE,
                "Treatment": MISSING_VALUE,
                "Replicate": replicates.get(treatment, MISSING_VALUE),
                "bamReads": f"{args.bam_dir}/{treatment}.sorted.rmdup.bam",
                "ControlID": control_id,
                "bamControl": bam_control,
                "Peaks": peak_path(treatment, args.peak_mode, args.peak_type),
                "PeakCaller": args.peak_type,
            }
        )

    return rows


def write_tsv(path: str | Path, rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", default="config/samples.yml", help="Path to samples.yml")
    parser.add_argument("--replicates", default="config/replicates.yml", help="Path to replicates.yml")
    parser.add_argument("--output", default="reports/chipqc/chipqc_sample_sheet.generated.tsv")
    parser.add_argument("--peak-mode", choices=["with_control", "without_control"], default="with_control")
    parser.add_argument("--peak-type", choices=["narrow", "broad"], default="narrow")
    parser.add_argument("--bam-dir", default="aligned_data")
    parser.add_argument("--no-control-value", default="NA")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = build_rows(args)
    write_tsv(args.output, rows)
    print(f"Wrote {len(rows)} rows to {Path(args.output)}")


if __name__ == "__main__":
    main()
