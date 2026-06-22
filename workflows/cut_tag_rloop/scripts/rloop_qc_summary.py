#!/usr/bin/env python3
"""Create R-loop/CUT&Tag-focused QC summaries from workflow outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_matrix(path: str):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        samples = [field for field in reader.fieldnames or [] if field not in {"peak_id", "chrom", "start", "end"}]
        rows = []
        for row in reader:
            rows.append([float(row[sample]) for sample in samples])
    return samples, np.asarray(rows, dtype=float)


def write_corr(path: str, samples: list[str], matrix: np.ndarray) -> np.ndarray:
    log_counts = np.log2(matrix + 1)
    corr = np.corrcoef(log_counts.T) if len(samples) > 1 and matrix.shape[0] > 1 else np.eye(len(samples))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *samples])
        for sample, values in zip(samples, corr):
            writer.writerow([sample, *[f"{value:.6f}" for value in values]])
    return corr


def plot_heatmap(path: str, samples: list[str], corr: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(max(4, len(samples) * 0.6), max(4, len(samples) * 0.6)))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(len(samples)), samples, rotation=90)
    ax.set_yticks(range(len(samples)), samples)
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_pca(path: str, samples: list[str], matrix: np.ndarray) -> None:
    log_counts = np.log2(matrix + 1)
    centered = log_counts.T - log_counts.T.mean(axis=0, keepdims=True)
    if centered.shape[0] >= 2 and centered.shape[1] >= 2:
        u, s, _ = np.linalg.svd(centered, full_matrices=False)
        coords = u[:, :2] * s[:2]
    else:
        coords = np.zeros((len(samples), 2))
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(coords[:, 0], coords[:, 1], s=36)
    for sample, x, y in zip(samples, coords[:, 0], coords[:, 1]):
        ax.text(x, y, sample, fontsize=7)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_peak_widths(path: str, peak_files: list[str]) -> None:
    widths = []
    for peak_file in peak_files:
        with open(peak_file) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    try:
                        widths.append(int(fields[2]) - int(fields[1]))
                    except ValueError:
                        pass
    fig, ax = plt.subplots(figsize=(5, 4))
    if widths:
        ax.hist(widths, bins=50, color="#3b6ea8")
    ax.set_xlabel("Peak width (bp)")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def concat_tables(paths: list[str], output: str, source_label: str = "source") -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    rows = []
    fields = [source_label]
    seen = {source_label}
    for path in paths:
        if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
            continue
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for field in reader.fieldnames or []:
                if field not in seen:
                    fields.append(field)
                    seen.add(field)
            for row in reader:
                row[source_label] = Path(path).name
                rows.append(row)

    with open(output, "w", newline="") as out:
        if not rows:
            out.write(f"{source_label}\tmessage\nnone\tNo input files available.\n")
            return
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str | None) -> float:
    if value is None or value == "" or value == "NA":
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")


def finite_median(values: list[float]) -> str:
    finite = [value for value in values if np.isfinite(value)]
    if not finite:
        return "NA"
    return f"{float(np.median(finite)):.8g}"


def ratio_table_summary(path: str) -> tuple[dict[str, str], list[dict[str, float | str]]]:
    rows = []
    ratios = []
    differences = []
    sensitive = 0
    sample = Path(path).name
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = str(row.get("name") or sample)
            treatment = as_float(row.get("treatment_signal"))
            rnaseh = as_float(row.get("rnaseh_signal"))
            ratio = as_float(row.get("signal_ratio"))
            diff = as_float(row.get("signal_difference"))
            ratios.append(ratio)
            differences.append(diff)
            classification = row.get("classification", "")
            if classification == "exploratory_rnaseh_sensitive_ratio":
                sensitive += 1
            rows.append(
                {
                    "source": Path(path).name,
                    "treatment_signal": treatment,
                    "rnaseh_signal": rnaseh,
                    "signal_ratio": ratio,
                    "classification": classification,
                }
            )
    total = len(rows)
    summary = {
        "source": Path(path).name,
        "source_type": "rnaseh_ratio",
        "total_regions": str(total),
        "sensitive_or_depleted_regions": str(sensitive),
        "fraction": f"{(sensitive / total) if total else 0:.6f}",
        "median_signal_ratio": finite_median(ratios),
        "median_signal_difference": finite_median(differences),
        "median_log2fc": "NA",
        "significant_fraction": "NA",
    }
    return summary, rows


def deseq2_table_summary(path: str) -> dict[str, str]:
    log2fc = []
    depleted = 0
    total = 0
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            log2fc.append(as_float(row.get("log2FoldChange")))
            if row.get("classification") == "rnaseh_depleted":
                depleted += 1
    fraction = (depleted / total) if total else 0
    return {
        "source": Path(path).name,
        "source_type": "rnaseh_deseq2",
        "total_regions": str(total),
        "sensitive_or_depleted_regions": str(depleted),
        "fraction": f"{fraction:.6f}",
        "median_signal_ratio": "NA",
        "median_signal_difference": "NA",
        "median_log2fc": finite_median(log2fc),
        "significant_fraction": f"{fraction:.6f}",
    }


def write_rnaseh_specificity_summary(path: str, ratio_tables: list[str], deseq2_tables: list[str]) -> list[dict[str, str]]:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    summaries = []
    for table in ratio_tables:
        if table and Path(table).exists() and Path(table).stat().st_size:
            summary, _ = ratio_table_summary(table)
            summaries.append(summary)
    for table in deseq2_tables:
        if table and Path(table).exists() and Path(table).stat().st_size:
            summaries.append(deseq2_table_summary(table))

    fields = [
        "source",
        "source_type",
        "total_regions",
        "sensitive_or_depleted_regions",
        "fraction",
        "median_signal_ratio",
        "median_signal_difference",
        "median_log2fc",
        "significant_fraction",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        if summaries:
            writer.writerows(summaries)
    return summaries


def plot_rnaseh_signal_scatter(path: str, ratio_tables: list[str]) -> None:
    points = []
    for table in ratio_tables:
        if table and Path(table).exists() and Path(table).stat().st_size:
            _, rows = ratio_table_summary(table)
            points.extend(rows)

    fig, ax = plt.subplots(figsize=(5, 5))
    if points:
        colors = [
            "#d64f3f" if row.get("classification") == "exploratory_rnaseh_sensitive_ratio" else "#4f78a8"
            for row in points
        ]
        x = [float(row["rnaseh_signal"]) for row in points]
        y = [float(row["treatment_signal"]) for row in points]
        ax.scatter(x, y, s=10, alpha=0.55, c=colors, edgecolors="none")
        finite_values = [value for value in [*x, *y] if np.isfinite(value)]
        if finite_values:
            low = min(0.0, min(finite_values))
            high = max(finite_values)
            ax.plot([low, high], [low, high], color="#444444", linewidth=1, linestyle="--")
            ax.set_xlim(left=low)
            ax.set_ylim(bottom=low)
    else:
        ax.text(0.5, 0.5, "No RNaseH ratio tables available", ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel("RNase H signal")
    ax.set_ylabel("No RNase H signal")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_rnaseh_depletion_fraction(path: str, summaries: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(max(5, len(summaries) * 0.7), 4))
    if summaries:
        labels = [row["source"].replace(".rnaseh_sensitive_ratio.tsv", "").replace(".rnaseh_depleted.deseq2.tsv", "") for row in summaries]
        values = [float(row["fraction"]) for row in summaries]
        colors = ["#d64f3f" if row["source_type"] == "rnaseh_ratio" else "#4f78a8" for row in summaries]
        ax.bar(range(len(values)), values, color=colors)
        ax.set_xticks(range(len(labels)), labels, rotation=90)
        ax.set_ylim(0, max(1.0, max(values) * 1.1))
    else:
        ax.text(0.5, 0.5, "No RNaseH depletion summaries available", ha="center", va="center", transform=ax.transAxes)
        ax.set_ylim(0, 1)
    ax.set_ylabel("RNaseH-sensitive/depleted fraction")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_blacklist_mito_summary(path: str, samples: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "metric", "value", "note"])
        for sample in samples:
            writer.writerow([sample, "blacklist_mito_filter", "applied_upstream_if_configured", "See alignment/filter logs for exact removed counts."])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--peaks", nargs="*", default=[])
    parser.add_argument("--frip", nargs="*", default=[])
    parser.add_argument("--spikein", default="")
    parser.add_argument("--rnaseh-ratio-tables", nargs="*", default=[])
    parser.add_argument("--rnaseh-ratio-summaries", nargs="*", default=[])
    parser.add_argument("--rnaseh-deseq2", nargs="*", default=[])
    parser.add_argument("--replicate-correlation", required=True)
    parser.add_argument("--replicate-correlation-heatmap", required=True)
    parser.add_argument("--sample-pca", required=True)
    parser.add_argument("--peak-width-distribution", required=True)
    parser.add_argument("--rnaseh-depletion-summary", required=True)
    parser.add_argument("--rnaseh-specificity-summary", required=True)
    parser.add_argument("--rnaseh-signal-scatter", required=True)
    parser.add_argument("--rnaseh-depletion-fraction", required=True)
    parser.add_argument("--spikein-summary", required=True)
    parser.add_argument("--blacklist-mito-summary", required=True)
    parser.add_argument("--frip-fragment-level", required=True)
    args = parser.parse_args()

    samples, matrix = read_matrix(args.counts)
    corr = write_corr(args.replicate_correlation, samples, matrix)
    plot_heatmap(args.replicate_correlation_heatmap, samples, corr)
    plot_pca(args.sample_pca, samples, matrix)
    plot_peak_widths(args.peak_width_distribution, args.peaks)
    concat_tables([*args.rnaseh_ratio_summaries, *args.rnaseh_deseq2], args.rnaseh_depletion_summary)
    specificity_rows = write_rnaseh_specificity_summary(
        args.rnaseh_specificity_summary,
        args.rnaseh_ratio_tables,
        args.rnaseh_deseq2,
    )
    plot_rnaseh_signal_scatter(args.rnaseh_signal_scatter, args.rnaseh_ratio_tables)
    plot_rnaseh_depletion_fraction(args.rnaseh_depletion_fraction, specificity_rows)
    concat_tables([args.spikein] if args.spikein else [], args.spikein_summary)
    concat_tables(args.frip, args.frip_fragment_level)
    write_blacklist_mito_summary(args.blacklist_mito_summary, samples)


if __name__ == "__main__":
    main()
