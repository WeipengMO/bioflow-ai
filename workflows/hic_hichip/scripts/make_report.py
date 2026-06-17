#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import html
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--qc", required=True)
    p.add_argument("--peak-sources", nargs="*", default=[])
    p.add_argument("--loops", nargs="*", default=[])
    p.add_argument("--outdir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--comparisons", default="")
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    samples = [s for s in args.samples.split(",") if s]
    comparisons = [c for c in args.comparisons.split(",") if c]
    qc_html = "<p>QC summary not available.</p>"
    qc = pd.DataFrame()
    if Path(args.qc).exists():
        qc = pd.read_csv(args.qc, sep="\t")
        qc_html = qc.to_html(index=False, escape=True)
    peak_html = "<p>Peak source table not available.</p>"
    if args.peak_sources and args.peak_sources[0] and Path(args.peak_sources[0]).exists():
        peak_html = pd.read_csv(args.peak_sources[0], sep="\t").to_html(index=False, escape=True)
    loop_counts = []
    for path in args.loops:
        if not path or not Path(path).exists():
            continue
        sample = Path(path).parent.name
        try:
            df = pd.read_csv(path, sep="\t")
            loop_counts.append({"sample": sample, "loop_count": len(df)})
        except Exception:
            loop_counts.append({"sample": sample, "loop_count": "NA"})
    loop_html = pd.DataFrame(loop_counts).to_html(index=False, escape=True) if loop_counts else "<p>No loop table available.</p>"
    warnings = [
        "HiChIP auto peaks are HiChIP-derived pseudo peaks and should be treated as fallback peaks, not independent matched ChIP-seq evidence.",
        "Hi-C loops are not mark-specific and should not be interpreted as enhancer-promoter loops without external annotation.",
        "Differential loop results without biological replicates are descriptive fold-change only, not formal statistical significance.",
    ]
    body = f"""
    <html><head><meta charset='utf-8'><title>BioFlowAI Hi-C / HiChIP report</title></head>
    <body>
    <h1>BioFlowAI Hi-C / HiChIP report</h1>
    <h2>Samples</h2><p>{html.escape(', '.join(samples))}</p>
    <h2>Comparisons</h2><p>{html.escape(', '.join(comparisons) or 'None')}</p>
    <h2>HiC-Pro QC summary</h2>{qc_html}
    <h2>Peak sources</h2>{peak_html}
    <h2>Loop counts</h2>{loop_html}
    <h2>Warnings</h2><ul>{''.join(f'<li>{html.escape(w)}</li>' for w in warnings)}</ul>
    <h2>Key output directories</h2>
    <ul>
      <li>Contact matrices: {html.escape(args.outdir)}/cool</li>
      <li>Loops: {html.escape(args.outdir)}/loops</li>
      <li>Differential loops: {html.escape(args.outdir)}/diffloop</li>
      <li>QC: {html.escape(args.outdir)}/qc</li>
    </ul>
    </body></html>
    """
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body)


if __name__ == "__main__":
    main()
