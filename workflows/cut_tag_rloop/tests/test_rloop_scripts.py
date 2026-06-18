#!/usr/bin/env python3
"""Lightweight script-level tests for CUT&Tag R-loop helpers."""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rnaseh_ratio_uses_pseudocount():
    rnaseh_signal = load_script("rnaseh_signal.py")
    ratio = rnaseh_signal.ratio_with_pseudocount(1.0, 0.0, 0.1)
    assert round(ratio, 6) == 11.0


def test_spikein_proper_pair_fragment_count_is_half_of_reads(monkeypatch):
    metrics = load_script("normalization_metrics.py")
    monkeypatch.setattr(metrics, "samtools_count", lambda args: 8)
    assert metrics.proper_pair_fragments("dummy.bam") == 4


def test_peak_universe_consensus_support():
    count_matrix = load_script("peak_count_matrix.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        peak1 = tmp / "rep1.bed"
        peak2 = tmp / "rep2.bed"
        out = tmp / "universe.bed"
        peak1.write_text("chr1\t10\t80\nchr1\t200\t280\n")
        peak2.write_text("chr1\t50\t120\nchr1\t500\t580\n")
        count_matrix.build_universe([str(peak1), str(peak2)], str(out), 50, 0, "consensus", 2)
        rows = out.read_text().strip().splitlines()
        assert rows == ["chr1\t10\t120\tpeak_1\t2\t."]


def test_parse_featurecounts_matrix():
    count_matrix = load_script("peak_count_matrix.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        fc = tmp / "featureCounts.txt"
        fc.write_text(
            "# Program:featureCounts\n"
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\tb.bam\n"
            "peak_1\tchr1\t11\t120\t.\t110\t5\t7\n"
            "peak_2\tchr2\t21\t70\t.\t50\t0\t3\n"
        )
        rows = [["chr1", "10", "120", "peak_1"], ["chr2", "20", "70", "peak_2"]]
        assert count_matrix.parse_featurecounts(str(fc), ["a", "b"], rows) == [[5, 7], [0, 3]]
