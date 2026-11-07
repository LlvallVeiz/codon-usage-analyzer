#!/usr/bin/env python3
"""Entry point: codon-usage analysis for a FASTA file.

Usage:
    python3 analyze.py path/to/sequences.fasta [--outdir out/]

Writes codon_table.txt, codon_table.csv and rscu_chart.txt into the
output directory and prints per-sequence + global summaries
(RSCU inputs, GC/GC1/GC2/GC3, Wright's Nc). Pure standard library.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from codon_usage import cli  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(cli.main(sys.argv[1:]))
