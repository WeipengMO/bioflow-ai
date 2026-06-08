#!/usr/bin/env python3
"""Render a project-specific HiC-Pro config from BioFlowAI YAML."""
from __future__ import annotations
import argparse
from pathlib import Path
import yaml


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--template", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text()) or {}
    genome = cfg.get("genome", {}) or {}
    hicpro = cfg.get("hicpro", {}) or {}
    outdir = Path(cfg.get("outdir", "results"))
    bowtie_prefix = Path(str(genome.get("bowtie2_index", "")))
    reference_genome = bowtie_prefix.name
    bowtie2_index_path = str(bowtie_prefix.parent) if str(bowtie_prefix.parent) != "." else str(bowtie_prefix)
    resolutions = hicpro.get("resolutions", [10000])
    bin_size = " ".join(str(int(x)) for x in resolutions)
    threads = (cfg.get("threads", {}) or {}).get("hicpro", 32)
    values = {
        "tmp_dir": str(outdir / "tmp" / "hicpro"),
        "logs_dir": str(outdir / "logs" / "hicpro"),
        "bowtie2_index_path": bowtie2_index_path,
        "reference_genome": reference_genome,
        "chrom_sizes": str(genome.get("chrom_sizes", "")),
        "restriction_fragments": str(genome.get("restriction_fragments", "")),
        "ligation_site": str(hicpro.get("ligation_site", "GATCGATC")),
        "min_mapq": str(int(hicpro.get("min_mapq", 30))),
        "pair1_ext": "_R1.fastq.gz",
        "pair2_ext": "_R2.fastq.gz",
        "bin_size": bin_size,
        "n_cpu": str(int(threads)),
        "rm_dup": "1" if bool(hicpro.get("remove_duplicates", True)) else "0",
    }
    template = Path(args.template).read_text()
    rendered = template
    for key, value in values.items():
        rendered = rendered.replace("{{" + key + "}}", value)
    missing = [part.split("}}", 1)[0] for part in rendered.split("{{")[1:]]
    if missing:
        raise SystemExit("Unresolved template placeholders: " + ", ".join(missing))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(rendered)


if __name__ == "__main__":
    main()
