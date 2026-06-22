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


def test_absolute_spikein_scale_factor_uses_constant_over_count():
    metrics = load_script("normalization_metrics.py")
    assert metrics.absolute_spikein_scale_factor(2500) == 4.0
    assert metrics.absolute_spikein_scale_factor(0) is None


def test_spikein_group_manifest_uses_anchor_column():
    metrics = load_script("normalization_metrics.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        manifest = Path(tmpdir) / "spikein_groups.tsv"
        manifest.write_text(
            "sample\tspikein_group\tspikein_anchor_sample\n"
            "target_sample\tspikein_anchor\tanchor_sample\n"
        )
        assert metrics.read_spikein_groups(str(manifest)) == {
            "target_sample": {
                "spikein_group": "spikein_anchor",
                "spikein_anchor_sample": "anchor_sample",
            }
        }


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


def test_spikein_matrix_not_written_without_metrics():
    count_matrix = load_script("peak_count_matrix.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        out = tmp / "peak_counts.spikein_normalized.tsv"
        wrote = count_matrix.write_spikein_matrix_if_requested(
            str(out),
            "",
            [["chr1", "10", "20", "peak_1"]],
            ["sample_a"],
            [[10]],
        )
        assert wrote is False
        assert not out.exists()


def test_spikein_matrix_uses_scale_factors():
    count_matrix = load_script("peak_count_matrix.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        metrics = tmp / "spikein_summary.tsv"
        out = tmp / "peak_counts.spikein_normalized.tsv"
        metrics.write_text(
            "sample\tmatched_anchor_spikein_scale_factor\n"
            "sample_a\t2\n"
            "sample_b\t0.5\n"
        )
        wrote = count_matrix.write_spikein_matrix_if_requested(
            str(out),
            str(metrics),
            [["chr1", "10", "20", "peak_1"]],
            ["sample_a", "sample_b"],
            [[10, 8]],
        )
        assert wrote is True
        assert out.read_text().splitlines() == [
            "peak_id\tchrom\tstart\tend\tsample_a\tsample_b",
            "peak_1\tchr1\t10\t20\t20\t4",
        ]


def test_frip_featurecounts_parser_sums_fragment_counts():
    frip = load_script("frip_score.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        fc = Path(tmpdir) / "featureCounts.txt"
        fc.write_text(
            "# Program:featureCounts\n"
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\n"
            "peak_1\tchr1\t11\t20\t.\t10\t3\n"
            "peak_2\tchr1\t31\t40\t.\t10\t5\n"
        )
        assert frip.parse_featurecounts_assigned(str(fc)) == 8


def test_rloop_qc_specificity_outputs():
    qc = load_script("rloop_qc_summary.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        ratio = tmp / "sample.rnaseh_sensitive_ratio.tsv"
        deseq = tmp / "group.rnaseh_depleted.deseq2.tsv"
        summary = tmp / "rnaseh_specificity_summary.tsv"
        scatter = tmp / "rnaseh_signal_scatter.pdf"
        fraction = tmp / "rnaseh_depletion_fraction.pdf"
        ratio.write_text(
            "chrom\tstart\tend\tname\ttreatment_signal\trnaseh_signal\tsignal_ratio\tsignal_difference\tpseudocount\tis_exploratory\tclassification\tscale_method\n"
            "chr1\t10\t20\tpeak_1\t4\t1\t3.727\t3\t0.1\ttrue\texploratory_rnaseh_sensitive_ratio\tspikein\n"
            "chr1\t30\t40\tpeak_2\t2\t2\t1\t0\t0.1\ttrue\tnot_sensitive\tspikein\n"
        )
        deseq.write_text(
            "peak_id\tchrom\tstart\tend\tlog2FoldChange\tpadj\tclassification\n"
            "peak_1\tchr1\t10\t20\t2\t0.01\trnaseh_depleted\n"
            "peak_2\tchr1\t30\t40\t0.2\t0.9\tnot_significant\n"
        )
        rows = qc.write_rnaseh_specificity_summary(str(summary), [str(ratio)], [str(deseq)])
        qc.plot_rnaseh_signal_scatter(str(scatter), [str(ratio)])
        qc.plot_rnaseh_depletion_fraction(str(fraction), rows)

        lines = summary.read_text().splitlines()
        assert lines[0].startswith("source\tsource_type\ttotal_regions")
        assert "sample.rnaseh_sensitive_ratio.tsv\trnaseh_ratio\t2\t1\t0.500000" in lines[1]
        assert scatter.stat().st_size > 0
        assert fraction.stat().st_size > 0
