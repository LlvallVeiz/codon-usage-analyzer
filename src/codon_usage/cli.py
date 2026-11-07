"""CLI: analyse a FASTA file and write codon-usage reports.

Outputs (into --outdir, default: current directory):
    codon_table.txt   fixed-width table (count, /1000, fraction, RSCU)
    codon_table.csv   the same table as CSV
    rscu_chart.txt    text bar chart of RSCU values
and a per-sequence + global summary printed to stdout.
"""

from __future__ import annotations

import argparse
import os
import sys

from . import analysis, fasta, report


def _summary_block(name: str, bundle: dict) -> str:
    gc = bundle["gc"]
    nc = bundle["nc"]
    lines = [
        f"--- {name} ---",
        f"sense codons : {bundle['sense_codons']}",
        f"stop codons  : {bundle['stop_codons']}"
        f"   skipped (ambiguous) : {bundle['skipped_codons']}"
        f"   leftover bases : {bundle['leftover_bases']}",
        f"GC           : {gc['gc']:.3f}   GC1 {gc['gc1']:.3f}  "
        f"GC2 {gc['gc2']:.3f}  GC3 {gc['gc3']:.3f}",
        f"Nc (Wright)  : {nc:.2f}" if nc is not None
        else "Nc (Wright)  : n/a (no usable family)",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="analyze", description="Codon usage + RSCU + GC3 + Nc "
                                    "for the CDS sequences in a FASTA")
    p.add_argument("fasta", help="input FASTA file")
    p.add_argument("--outdir", default=".",
                   help="where to write the tables (default: .)")
    args = p.parse_args(argv)

    records = fasta.read_fasta(args.fasta)
    if not records:
        print("no records found", file=sys.stderr)
        return 1

    per_seq = []
    for header, seq in records:
        bundle = analysis.analyse_sequence(seq)
        bundle["name"] = fasta.record_id(header)
        per_seq.append(bundle)
        print(_summary_block(bundle["name"], bundle))

    pooled_counts = analysis.merge_counts(b["counts"] for b in per_seq)
    pooled = {
        "counts": pooled_counts,
        "skipped_codons": sum(b["skipped_codons"] for b in per_seq),
        "leftover_bases": sum(b["leftover_bases"] for b in per_seq),
        "stop_codons": sum(b["stop_codons"] for b in per_seq),
        "sense_codons": analysis.sense_total(pooled_counts),
        "rscu": analysis.rscu(pooled_counts),
        "fraction": analysis.family_fraction(pooled_counts),
        "per_thousand": analysis.per_thousand(pooled_counts),
        "gc": analysis.gc_by_position(pooled_counts),
        "nc": analysis.effective_number_of_codons(pooled_counts),
    }
    print(_summary_block(f"GLOBAL (pooled, {len(records)} sequences)",
                         pooled))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "codon_table.txt"), "w") as fh:
        fh.write(report.table_text(pooled_counts, pooled["fraction"],
                                   pooled["per_thousand"],
                                   pooled["rscu"]))
    with open(os.path.join(args.outdir, "codon_table.csv"), "w") as fh:
        fh.write(report.table_csv(pooled_counts, pooled["fraction"],
                                  pooled["per_thousand"],
                                  pooled["rscu"]))
    with open(os.path.join(args.outdir, "rscu_chart.txt"), "w") as fh:
        fh.write(report.rscu_chart(pooled["rscu"]))
    print(f"tables written to {args.outdir}/ "
          f"(codon_table.txt, codon_table.csv, rscu_chart.txt)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
