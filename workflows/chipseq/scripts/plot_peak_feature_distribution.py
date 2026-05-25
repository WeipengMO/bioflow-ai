from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "fontconfig-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MAJOR_FEATURE_ORDER = ["promoter", "exon", "intron", "intergenic", "TTS", "other"]
FEATURE_COLORS = {
    "promoter": "#D55E00",
    "exon": "#E69F00",
    "intron": "#56B4E9",
    "intergenic": "#009E73",
    "TTS": "#CC79A7",
    "other": "#7F7F7F",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--summary-inputs", nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary-tsv", required=True)
    args = parser.parse_args()
    if bool(args.input) == bool(args.summary_inputs):
        parser.error("Specify exactly one of --input or --summary-inputs.")
    return args


def normalize_annotation(annotation: str) -> str:
    token = annotation.strip()
    if not token:
        return "other"
    token = token.split(" ", 1)[0]
    token = token.split("(", 1)[0]
    token = token.strip()
    token_lower = token.lower()
    if token_lower.startswith("promoter"):
        return "promoter"
    if token_lower.startswith("exon"):
        return "exon"
    if token_lower.startswith("intron"):
        return "intron"
    if token_lower.startswith("intergenic"):
        return "intergenic"
    if token == "TTS" or token_lower.startswith("tts"):
        return "TTS"
    return "other"


def load_counts(path: str) -> Counter:
    counts: Counter = Counter()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "Annotation" not in (reader.fieldnames or []):
            raise ValueError(f"Missing required 'Annotation' column in {path}")
        for row in reader:
            counts[normalize_annotation(row.get("Annotation", ""))] += 1
    return counts


def ordered_counts(counts: Counter) -> tuple[list[str], list[int]]:
    return MAJOR_FEATURE_ORDER.copy(), [counts.get(feature, 0) for feature in MAJOR_FEATURE_ORDER]


def render_pie(labels: list[str], values: list[int], output_path: str) -> None:
    total = sum(values)
    fig, ax = plt.subplots(figsize=(7, 7), dpi=160)
    if total == 0:
        ax.text(0.5, 0.5, "No annotated peaks", ha="center", va="center", fontsize=13)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        return

    nonzero = [(label, value) for label, value in zip(labels, values) if value > 0]
    plot_labels = [label for label, _ in nonzero]
    plot_values = [value for _, value in nonzero]
    colors = [FEATURE_COLORS[label] for label in plot_labels]
    display_labels = [f"{label} ({value:,})" for label, value in nonzero]

    ax.pie(
        plot_values,
        labels=display_labels,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 2 else "",
        startangle=90,
        counterclock=False,
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
        textprops={"fontsize": 10},
    )
    ax.set_title(f"Peak feature distribution (n={total:,})", fontsize=13)
    ax.axis("equal")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def write_summary(labels: list[str], values: list[int], output_path: str) -> None:
    total = sum(values)
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature", "count", "fraction", "percent"])
        for label, value in zip(labels, values):
            fraction = value / total if total else 0
            writer.writerow([label, value, f"{fraction:.6f}", f"{fraction * 100:.2f}"])


def load_summary(path: str) -> dict[str, int]:
    counts = {feature: 0 for feature in MAJOR_FEATURE_ORDER}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"feature", "count"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required column(s) in {path}: {', '.join(sorted(missing))}")
        for row in reader:
            feature = row["feature"]
            if feature not in counts:
                continue
            counts[feature] = int(row["count"])
    return counts


def target_from_summary_path(path: str) -> str:
    return Path(path).parent.name


def load_aggregate_rows(paths: list[str]) -> list[dict[str, object]]:
    rows = []
    seen_targets = set()
    for path in paths:
        target = target_from_summary_path(path)
        if target in seen_targets:
            raise ValueError(f"Duplicate HOMER target name in summary inputs: {target}")
        seen_targets.add(target)
        counts = load_summary(path)
        total = sum(counts.values())
        rows.append({"target": target, "total_peaks": total, **counts})
    return rows


def write_aggregate_summary(rows: list[dict[str, object]], output_path: str) -> None:
    header = (
        ["target", "total_peaks"]
        + [f"{feature}_count" for feature in MAJOR_FEATURE_ORDER]
        + [f"{feature}_fraction" for feature in MAJOR_FEATURE_ORDER]
    )
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            total = int(row["total_peaks"])
            record = {"target": row["target"], "total_peaks": total}
            for feature in MAJOR_FEATURE_ORDER:
                count = int(row[feature])
                record[f"{feature}_count"] = count
                record[f"{feature}_fraction"] = f"{(count / total if total else 0):.6f}"
            writer.writerow(record)


def render_stacked_bar(rows: list[dict[str, object]], output_path: str) -> None:
    row_count = max(len(rows), 1)
    fig_height = max(3.2, min(18.0, 1.2 + row_count * 0.42))
    fig, ax = plt.subplots(figsize=(10, fig_height), dpi=160)

    if not rows:
        ax.text(0.5, 0.5, "No peak feature summaries", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        return

    y_positions = list(range(len(rows)))
    left = [0.0] * len(rows)
    totals = [int(row["total_peaks"]) for row in rows]

    for feature in MAJOR_FEATURE_ORDER:
        values = [
            (int(row[feature]) / total * 100) if total else 0.0
            for row, total in zip(rows, totals)
        ]
        ax.barh(
            y_positions,
            values,
            left=left,
            color=FEATURE_COLORS[feature],
            edgecolor="white",
            linewidth=0.5,
            label=feature,
        )
        left = [base + value for base, value in zip(left, values)]

    labels = [f"{row['target']} (n={int(row['total_peaks']):,})" for row in rows]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Peaks (%)")
    ax.set_title("Peak feature distribution")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.summary_inputs:
        rows = load_aggregate_rows(args.summary_inputs)
        write_aggregate_summary(rows, args.summary_tsv)
        render_stacked_bar(rows, args.output)
        for row in rows:
            print(f"{row['target']}\t{row['total_peaks']}")
    else:
        counts = load_counts(args.input)
        labels, values = ordered_counts(counts)
        render_pie(labels, values, args.output)
        write_summary(labels, values, args.summary_tsv)
        for label, value in zip(labels, values):
            print(f"{label}\t{value}")


if __name__ == "__main__":
    main()
