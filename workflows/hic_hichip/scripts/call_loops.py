#!/usr/bin/env python3
"""Dispatch loop calling and standardize outputs to BEDPE/TSV/anchors.

This wrapper intentionally fails with actionable messages when an external caller
is not installed. It also supports precomputed loop import through
loop_calling.precomputed_loops in config.yml for testing or reanalysis.
"""
from __future__ import annotations
import argparse
import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import yaml
import pandas as pd

BEDPE_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def as_bool(value, default=False):
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "yes", "y", "1", "on"}:
        return True
    if normalized in {"false", "f", "no", "n", "0", "off", ""}:
        return False
    return bool(default)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--assay", required=True)
    p.add_argument("--target", default="")
    p.add_argument("--caller", required=True)
    p.add_argument("--resolution", required=True, type=int)
    p.add_argument("--fdr", required=True)
    p.add_argument("--valid-pairs", required=True)
    p.add_argument("--cool", required=True)
    p.add_argument("--mcool", required=True)
    p.add_argument("--peaks", default="")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--bedpe", required=True)
    p.add_argument("--tsv", required=True)
    p.add_argument("--anchors", required=True)
    return p.parse_args()


def run(cmd, cwd=None):
    print("[call_loops]", " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def normalize_loop_table(path: Path, sample: str, caller: str, qvalue_default="NA") -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"Loop caller did not produce a non-empty loop file: {path}")
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path, sep="\t", header=None)
    lower = {str(c).lower(): c for c in df.columns}
    def col(*names):
        for name in names:
            if name in lower:
                return lower[name]
        return None
    if set(BEDPE_COLUMNS).issubset(set(df.columns)):
        out = df.copy()
    else:
        c1 = col("chrom1", "chr1", "chr_a", "chromosome1")
        s1 = col("start1", "x1", "start_a", "bin1_start")
        e1 = col("end1", "x2", "end_a", "bin1_end")
        c2 = col("chrom2", "chr2", "chr_b", "chromosome2")
        s2 = col("start2", "y1", "start_b", "bin2_start")
        e2 = col("end2", "y2", "end_b", "bin2_end")
        if None in [c1, s1, e1, c2, s2, e2]:
            # common 6-column BEDPE without header
            if df.shape[1] >= 6:
                out = df.iloc[:, :6].copy()
                out.columns = BEDPE_COLUMNS
            else:
                raise ValueError(f"Cannot infer BEDPE columns from {path}")
        else:
            out = df[[c1, s1, e1, c2, s2, e2]].copy()
            out.columns = BEDPE_COLUMNS
    for c in ["start1", "end1", "start2", "end2"]:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    out = out.dropna(subset=BEDPE_COLUMNS)
    out["loop_id"] = [f"{sample}_loop_{i+1}" for i in range(len(out))]
    out["score"] = df[col("score", "count", "observed", "balanced.avg")].values[:len(out)] if col("score", "count", "observed", "balanced.avg") else "NA"
    qcol = col("qvalue", "q.value", "qval", "fdr", "p.adj", "padj")
    out["qvalue"] = df[qcol].values[:len(out)] if qcol else qvalue_default
    out["sample"] = sample
    out["caller"] = caller
    return out[[*BEDPE_COLUMNS, "loop_id", "score", "qvalue", "sample", "caller"]]


def write_standard_outputs(df: pd.DataFrame, bedpe: Path, tsv: Path, anchors: Path):
    bedpe.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv, sep="\t", index=False)
    df.to_csv(bedpe, sep="\t", index=False)
    a1 = df[["chrom1", "start1", "end1"]].copy(); a1.columns = ["chrom", "start", "end"]
    a2 = df[["chrom2", "start2", "end2"]].copy(); a2.columns = ["chrom", "start", "end"]
    anchors_df = pd.concat([a1, a2], ignore_index=True).drop_duplicates().sort_values(["chrom", "start", "end"])
    anchors_df.to_csv(anchors, sep="\t", header=False, index=False)


def import_precomputed(cfg, sample):
    pre = (cfg.get("loop_calling", {}) or {}).get("precomputed_loops", {}) or {}
    if isinstance(pre, str) and pre:
        return Path(pre)
    if isinstance(pre, dict) and sample in pre:
        return Path(pre[sample])
    return None


def call_mustache(args, loop_cfg, tmpdir):
    exe = loop_cfg.get("mustache_bin", "mustache")
    if shutil.which(exe) is None:
        raise SystemExit(f"mustache executable not found: {exe}. Install it or set loop_calling.precomputed_loops.")
    out = Path(tmpdir) / "mustache_loops.tsv"
    run([exe, "-f", args.cool, "-r", str(args.resolution), "-pt", str(args.fdr), "-p", str(args.threads), "-o", str(out)])
    return out


def call_cooltools(args, loop_cfg, tmpdir):
    exe = loop_cfg.get("cooltools_bin", "cooltools")
    if shutil.which(exe) is None:
        raise SystemExit(f"cooltools executable not found: {exe}. Install it or set loop_calling.precomputed_loops.")
    out = Path(tmpdir) / "cooltools_dots.tsv"
    # Users may need to tune expected/dots parameters for their genome/binning.
    expected = Path(tmpdir) / "expected.tsv"
    run([exe, "expected-cis", args.cool, "-o", str(expected), "--nproc", str(args.threads)])
    try:
        run([exe, "dots", args.cool, str(expected), "-o", str(out), "--nproc", str(args.threads)])
    except subprocess.CalledProcessError:
        if not as_bool(loop_cfg.get("allow_empty", True), True):
            raise
        print(
            "[call_loops] WARNING: cooltools dots failed. Writing an empty loop table because "
            "loop_calling.allow_empty is true. This is expected for coarse-resolution or low-depth smoke tests."
        )
        pd.DataFrame(columns=BEDPE_COLUMNS).to_csv(out, sep="\t", index=False)
    return out


def executable_path(path_or_name):
    if not path_or_name:
        return None
    candidate = Path(str(path_or_name)).expanduser()
    if candidate.exists():
        return str(candidate.resolve())
    found = shutil.which(str(path_or_name))
    return found


def resolve_fithichip_script(loop_cfg):
    workflow_dir = Path(__file__).resolve().parents[1]
    candidates = [
        loop_cfg.get("fithichip_script", ""),
        os.environ.get("FITHICHIP_SCRIPT", ""),
        Path(os.environ.get("FITHICHIP_HOME", "")) / "FitHiChIP_HiCPro.sh" if os.environ.get("FITHICHIP_HOME") else "",
        Path.cwd() / "external" / "FitHiChIP" / "FitHiChIP_HiCPro.sh",
        workflow_dir / "external" / "FitHiChIP" / "FitHiChIP_HiCPro.sh",
        "FitHiChIP_HiCPro.sh",
    ]
    for candidate in candidates:
        resolved = executable_path(candidate)
        if resolved:
            return resolved
    searched = [str(c) for c in candidates if str(c)]
    raise SystemExit(
        "FitHiChIP runner not found. Run scripts/install_fithichip.sh --prefix tools/FitHiChIP "
        "and set FITHICHIP_SCRIPT in config/run_snakemake.env, or set loop_calling.fithichip_script. "
        f"Searched: {', '.join(searched)}"
    )


def call_fithichip(args, loop_cfg, tmpdir):
    script = resolve_fithichip_script(loop_cfg)
    if not args.peaks:
        raise SystemExit(
            "FitHiChIP requires peaks for the default peak-aware mode. "
            "Provide peak_bed in samples.tsv or loop_calling.external_peaks in config.yml."
        )
    cfg_path = Path(tmpdir) / f"{args.sample}.fithichip.config"
    outdir = Path(tmpdir) / "fithichip_out"
    cfg_path.write_text("\n".join([
        f"ValidPairs={Path(args.valid_pairs).resolve()}",
        f"PeakFile={Path(args.peaks).resolve()}",
        f"OutDir={outdir.resolve()}",
        f"PREFIX={args.sample}",
        f"Resolution={args.resolution}",
        f"FDRThr={args.fdr}",
        "UseP2PBackgrnd=0",
        "BiasType=Coverage",
        "MergeInt=1",
        "" ]))
    run(["bash", script, "-C", str(cfg_path)])
    candidates = list(outdir.glob("**/*FitHiChIP*.interactions_FitHiC_Q*.bed")) + list(outdir.glob("**/*loops*.bedpe")) + list(outdir.glob("**/*.bedpe"))
    if not candidates:
        raise SystemExit(f"FitHiChIP finished but no loop BEDPE-like output was found under {outdir}")
    return sorted(candidates)[0]


def main():
    args = parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text()) or {}
    loop_cfg = cfg.get("loop_calling", {}) or {}
    precomputed = import_precomputed(cfg, args.sample)
    with tempfile.TemporaryDirectory(prefix=f"{args.sample}.loops.") as tmpdir:
        if precomputed:
            loop_file = precomputed
        elif args.caller == "mustache":
            loop_file = call_mustache(args, loop_cfg, tmpdir)
        elif args.caller == "cooltools":
            loop_file = call_cooltools(args, loop_cfg, tmpdir)
        elif args.caller == "fithichip":
            loop_file = call_fithichip(args, loop_cfg, tmpdir)
        else:
            raise SystemExit(f"Unsupported loop caller: {args.caller}")
        df = normalize_loop_table(Path(loop_file), args.sample, args.caller, qvalue_default=args.fdr)
        write_standard_outputs(df, Path(args.bedpe), Path(args.tsv), Path(args.anchors))


if __name__ == "__main__":
    main()
