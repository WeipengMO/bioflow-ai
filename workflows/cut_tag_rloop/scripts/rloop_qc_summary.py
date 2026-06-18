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
    parser.add_argument("--rnaseh-ratio-summaries", nargs="*", default=[])
    parser.add_argument("--rnaseh-deseq2", nargs="*", default=[])
    parser.add_argument("--replicate-correlation", required=True)
    parser.add_argument("--replicate-correlation-heatmap", required=True)
    parser.add_argument("--sample-pca", required=True)
    parser.add_argument("--peak-width-distribution", required=True)
    parser.add_argument("--rnaseh-depletion-summary", required=True)
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
    concat_tables([args.spikein] if args.spikein else [], args.spikein_summary)
    concat_tables(args.frip, args.frip_fragment_level)
    write_blacklist_mito_summary(args.blacklist_mito_summary, samples)


if __name__ == "__main__":
    main()
