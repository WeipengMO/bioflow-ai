#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import html
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--qc", required=True)
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
    if Path(args.qc).exists():
        qc_html = pd.read_csv(args.qc, sep="\t").to_html(index=False, escape=True)
    body = f"""
    <html><head><meta charset='utf-8'><title>BioFlowAI Hi-C / HiChIP report</title></head>
    <body>
    <h1>BioFlowAI Hi-C / HiChIP report</h1>
    <h2>Samples</h2><p>{html.escape(', '.join(samples))}</p>
    <h2>Comparisons</h2><p>{html.escape(', '.join(comparisons) or 'None')}</p>
    <h2>HiC-Pro QC summary</h2>{qc_html}
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
