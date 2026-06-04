#!/usr/bin/env python3
"""Filter candidate lncRNAs by ORF length and optional coding-potential results."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


STOP_CODONS = {"TAA", "TAG", "TGA"}


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_attrs(attr: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr)}


def transcript_id(fields: List[str]) -> str | None:
    return parse_attrs(fields[8]).get("transcript_id")


def read_fasta(path: str) -> Dict[str, str]:
    seqs = {}
    name = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper().replace("U", "T")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper().replace("U", "T")
    return seqs


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def longest_orf_nt_single_strand(seq: str) -> int:
    longest = 0
    for frame in range(3):
        starts = []
        for pos in range(frame, len(seq) - 2, 3):
            codon = seq[pos : pos + 3]
            if codon == "ATG":
                starts.append(pos)
            if codon in STOP_CODONS and starts:
                stop_end = pos + 3
                longest = max(longest, max(stop_end - start for start in starts))
                starts = []
        if starts:
            longest = max(longest, max(len(seq) - start for start in starts))
    return longest


def longest_orf_nt(seq: str) -> int:
    return max(longest_orf_nt_single_strand(seq), longest_orf_nt_single_strand(reverse_complement(seq)))


def read_gtf_by_tx(path: str) -> Dict[str, List[str]]:
    records = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            tx_id = transcript_id(fields)
            if tx_id:
                records[tx_id].append(line)
    return records


def load_bad_ids(path: str | None) -> set[str]:
    if not path:
        return set()
    bad = set()
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            bad.add(line.split()[0])
    return bad


def main() -> None:
    parser = argparse.ArgumentParser(description="Create high-confidence lncRNA files from candidate lncRNAs.")
    parser.add_argument("--candidate-gtf", required=True)
    parser.add_argument("--candidate-fasta", required=True)
    parser.add_argument("--candidate-summary", required=True)
    parser.add_argument("--output-gtf", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--orf-metrics", required=True)
    parser.add_argument("--max-orf-aa", type=int, default=100)
    parser.add_argument("--exclude-transcripts", default=None, help="Optional newline-delimited transcript IDs to remove.")
    args = parser.parse_args()

    seqs = read_fasta(args.candidate_fasta)
    records_by_tx = read_gtf_by_tx(args.candidate_gtf)
    bad_ids = load_bad_ids(args.exclude_transcripts)

    orf_rows = []
    keep = set()
    for tx_id, seq in sorted(seqs.items()):
        orf_nt = longest_orf_nt(seq)
        orf_aa = orf_nt // 3
        coding_by_orf = orf_aa > args.max_orf_aa
        keep_tx = tx_id not in bad_ids and not coding_by_orf
        if keep_tx:
            keep.add(tx_id)
        orf_rows.append((tx_id, len(seq), orf_nt, orf_aa, int(coding_by_orf), int(keep_tx)))

    Path(args.output_gtf).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_gtf, "w") as out:
        for tx_id in sorted(keep):
            for line in records_by_tx.get(tx_id, []):
                out.write(line)

    with open(args.output_fasta, "w") as out:
        for tx_id in sorted(keep):
            seq = seqs[tx_id]
            out.write(f">{tx_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")

    with open_text(args.candidate_summary) as handle, open(args.output_summary, "w") as out:
        header = handle.readline()
        out.write(header.rstrip("\n") + "\tmax_orf_aa\n")
        for line in handle:
            if not line.strip():
                continue
            tx_id = line.split("\t", 1)[0]
            if tx_id in keep:
                orf_aa = next(row[3] for row in orf_rows if row[0] == tx_id)
                out.write(line.rstrip("\n") + f"\t{orf_aa}\n")

    with open(args.orf_metrics, "w") as out:
        out.write("transcript_id\ttranscript_length\tmax_orf_nt\tmax_orf_aa\tcoding_by_orf\tkept\n")
        for row in orf_rows:
            out.write("\t".join(map(str, row)) + "\n")

    if not keep:
        raise SystemExit("No high-confidence lncRNAs remained after coding-potential filtering.")


if __name__ == "__main__":
    main()
