"""Codon-usage statistics: counts, RSCU, GC/GC3, and Wright's Nc.

Definitions (see docs/RESEARCH.md for sources and discussion):

* **Codon counts** — complete codons in frame 0. Codons containing
  non-ACGT characters are skipped and reported as `skipped_codons`;
  trailing partial bases are reported as `leftover_bases`.
* **Per thousand** — count / total sense codons × 1000.
* **Fraction** — a codon's share within its synonymous family.
* **RSCU** (Sharp & Li 1987) — observed count divided by the count
  expected if all synonyms of the amino acid were used equally:
  RSCU_ij = n_i · x_ij / Σ_j x_ij. RSCU = 1 means no preference;
  it is defined for the 61 sense codons (single-codon families Met
  and Trp trivially give 1.0). Stop codons are counted but excluded
  from RSCU, GC3 and Nc.
* **GC / GC1 / GC2 / GC3** — (G+C) fraction overall and at each codon
  position; position statistics use sense codons only.
* **Nc** (Wright 1990) — the effective number of codons:
  Nc = 2 + 9/F̄2 + 1/F̄3 + 5/F̄4 + 3/F̄6, where F̄k is the mean
  homozygosity of the k-fold degenerate families and each family's
  F = (n·Σp² − 1)/(n − 1) with Wright's sampling-without-replacement
  correction. Nc ranges from 20 (one codon per amino acid) to 61
  (uniform usage). Families with fewer than 2 observations cannot
  yield F and are excluded from their class mean; if an entire class
  is unobserved, its F̄ is replaced by the mean of the available
  class means — a documented simplification (see RESEARCH.md).
"""

from __future__ import annotations

from .genetic_code import CODON_TABLE, DEGENERACY_CLASSES, FAMILIES


def count_codons(sequence: str) -> dict:
    """Count complete codons of one sequence (frame 0).

    Returns {"counts": {codon: n}, "skipped_codons": int,
             "leftover_bases": int, "stop_codons": int}.
    """
    sequence = sequence.upper()
    counts: dict[str, int] = {}
    skipped = 0
    usable = len(sequence) - (len(sequence) % 3)
    for i in range(0, usable, 3):
        codon = sequence[i:i + 3]
        if codon in CODON_TABLE:
            counts[codon] = counts.get(codon, 0) + 1
        else:
            skipped += 1
    stops = sum(counts.get(c, 0) for c in ("TAA", "TAG", "TGA"))
    return {"counts": counts, "skipped_codons": skipped,
            "leftover_bases": len(sequence) % 3, "stop_codons": stops}


def merge_counts(list_of_counts: list[dict[str, int]]) -> dict[str, int]:
    merged: dict[str, int] = {}
    for counts in list_of_counts:
        for codon, n in counts.items():
            merged[codon] = merged.get(codon, 0) + n
    return merged


def sense_total(counts: dict[str, int]) -> int:
    """Total sense (non-stop) codons."""
    return sum(n for c, n in counts.items() if CODON_TABLE[c] != "*")


def rscu(counts: dict[str, int]) -> dict[str, float]:
    """RSCU for every sense codon (0.0 for unobserved families)."""
    out: dict[str, float] = {}
    for aa, codons in FAMILIES.items():
        family_total = sum(counts.get(c, 0) for c in codons)
        k = len(codons)
        for codon in codons:
            if family_total == 0:
                out[codon] = 0.0
            else:
                out[codon] = k * counts.get(codon, 0) / family_total
    return out


def family_fraction(counts: dict[str, int]) -> dict[str, float]:
    """Each sense codon's share within its synonymous family."""
    out: dict[str, float] = {}
    for aa, codons in FAMILIES.items():
        family_total = sum(counts.get(c, 0) for c in codons)
        for codon in codons:
            out[codon] = (counts.get(codon, 0) / family_total
                          if family_total else 0.0)
    return out


def per_thousand(counts: dict[str, int]) -> dict[str, float]:
    """Codon frequency per 1000 sense codons (stops: per 1000 stops-free
    total as well — they are reported against the same denominator so
    the column sums to 1000 over sense codons)."""
    total = sense_total(counts)
    if total == 0:
        return {c: 0.0 for c in CODON_TABLE}
    return {c: counts.get(c, 0) / total * 1000.0 for c in CODON_TABLE}


def gc_by_position(counts: dict[str, int]) -> dict[str, float]:
    """GC fraction overall and per codon position (sense codons only
    for the positional values; `gc` uses all counted codons)."""
    pos = [{}, {}, {}]  # per position: {"GC": n, "AT": n}
    total_gc = total_bases = 0
    for codon, n in counts.items():
        for base in codon:
            total_bases += n
            if base in "GC":
                total_gc += n
        if CODON_TABLE[codon] == "*":
            continue
        for i, base in enumerate(codon):
            key = "GC" if base in "GC" else "AT"
            pos[i][key] = pos[i].get(key, 0) + n
    result = {"gc": (total_gc / total_bases) if total_bases else 0.0}
    for i in range(3):
        gc, at = pos[i].get("GC", 0), pos[i].get("AT", 0)
        result[f"gc{i + 1}"] = gc / (gc + at) if (gc + at) else 0.0
    return result


def _family_homozygosity(codon_counts: list[int]) -> float | None:
    """Wright's F for one amino-acid family; None if n < 2.

    F is deliberately NOT floored at 1/k: with Wright's correction,
    sampling noise can legitimately push the estimate below the
    uniform-usage value for small families, and flooring it would
    bias Nc upward for short genes. The [20, 61] clamp on Nc itself
    absorbs the extremes. (F <= 1 holds automatically since
    sum(p^2) <= 1; the min() is a numerical safeguard.)
    """
    n = sum(codon_counts)
    if n < 2:
        return None
    sum_p2 = sum((c / n) ** 2 for c in codon_counts)
    f = (n * sum_p2 - 1.0) / (n - 1.0)
    return min(1.0, f)


def effective_number_of_codons(counts: dict[str, int]) -> float | None:
    """Wright's Nc from a codon-count dict; None if no class is usable."""
    class_means: dict[int, float] = {}
    for degeneracy, aas in sorted(DEGENERACY_CLASSES.items()):
        if degeneracy == 1:
            continue  # Met/Trp enter via the constant term (2)
        fs = []
        for aa in aas:
            f = _family_homozygosity(
                [counts.get(c, 0) for c in FAMILIES[aa]])
            if f is not None:
                fs.append(f)
        if fs:
            class_means[degeneracy] = sum(fs) / len(fs)
    if not class_means:
        return None
    fallback = sum(class_means.values()) / len(class_means)

    def term(degeneracy: int, n_families: int) -> float:
        f = class_means.get(degeneracy, fallback)
        if f <= 0.0:
            # F = 0 arises when every observed codon in the class was
            # seen at most once: no bias is measurable, so the class
            # contributes its maximum (n_families x degeneracy).
            return float(n_families * degeneracy)
        return n_families / f

    nc = 2.0 + term(2, 9) + term(3, 1) + term(4, 5) + term(6, 3)
    return min(61.0, max(20.0, nc))


def analyse_sequence(sequence: str) -> dict:
    """Full per-sequence analysis bundle."""
    counted = count_codons(sequence)
    counts = counted["counts"]
    return {
        "counts": counts,
        "skipped_codons": counted["skipped_codons"],
        "leftover_bases": counted["leftover_bases"],
        "stop_codons": counted["stop_codons"],
        "sense_codons": sense_total(counts),
        "rscu": rscu(counts),
        "fraction": family_fraction(counts),
        "per_thousand": per_thousand(counts),
        "gc": gc_by_position(counts),
        "nc": effective_number_of_codons(counts),
    }
