"""Report rendering: codon table (text + CSV) and an RSCU bar chart."""

from __future__ import annotations

import csv
import io

from .genetic_code import AMINO_ACID_NAMES, CODON_TABLE, FAMILIES

_HEADER = (f"{'Codon':<6}{'AA':<10}{'Count':>7}{'/1000':>9}"
           f"{'Fraction':>10}{'RSCU':>8}")


def table_rows(counts, frac, per1000, rscu_values):
    """Yield (codon, aa_label, count, per1000, fraction, rscu|None)."""
    for aa in sorted(FAMILIES, key=lambda a: AMINO_ACID_NAMES[a]):
        for codon in FAMILIES[aa]:
            yield (codon, f"{aa} ({AMINO_ACID_NAMES[aa]})",
                   counts.get(codon, 0), per1000[codon],
                   frac[codon], rscu_values[codon])
    for codon in ("TAA", "TAG", "TGA"):
        yield (codon, "* (Stop)", counts.get(codon, 0),
               per1000[codon], 0.0, None)


def table_text(counts, frac, per1000, rscu_values) -> str:
    lines = [_HEADER, "-" * len(_HEADER)]
    current_aa = None
    for codon, aa_label, n, pt, fr, rs in table_rows(
            counts, frac, per1000, rscu_values):
        if aa_label != current_aa:
            lines.append("")
            current_aa = aa_label
        rscu_s = f"{rs:8.2f}" if rs is not None else "       -"
        lines.append(f"{codon:<6}{aa_label:<10}{n:>7}{pt:>9.1f}"
                     f"{fr:>10.3f}{rscu_s}")
    return "\n".join(lines) + "\n"


def table_csv(counts, frac, per1000, rscu_values) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["codon", "amino_acid", "count", "per_thousand",
                     "fraction_within_family", "rscu"])
    for codon, aa_label, n, pt, fr, rs in table_rows(
            counts, frac, per1000, rscu_values):
        aa = CODON_TABLE[codon]
        writer.writerow([codon, aa, n, f"{pt:.2f}", f"{fr:.4f}",
                         "" if rs is None else f"{rs:.4f}"])
    return buf.getvalue()


def rscu_chart(rscu_values, width: int = 40) -> str:
    """Text bar chart of RSCU values, grouped by amino acid.

    Bars are scaled so RSCU = family degeneracy (a codon taking the
    whole family) reaches `width`; a marker at RSCU 1.0 would sit at
    width/degeneracy — instead each line is labelled with its value,
    which keeps the chart honest across families of different sizes.
    """
    lines = ["RSCU by codon (1.00 = no preference within the family)",
             ""]
    for aa in sorted(FAMILIES, key=lambda a: AMINO_ACID_NAMES[a]):
        codons = FAMILIES[aa]
        lines.append(f"{AMINO_ACID_NAMES[aa]} ({aa}):")
        for codon in codons:
            value = rscu_values[codon]
            bar_len = int(round(value / len(codons) * width))
            lines.append(f"  {codon}  {'#' * bar_len:<{width}} {value:5.2f}")
        lines.append("")
    return "\n".join(lines)
