"""Tests with hand-computed expectations for every statistic."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from codon_usage import analysis  # noqa: E402
from codon_usage.genetic_code import (  # noqa: E402
    CODON_TABLE, DEGENERACY_CLASSES, FAMILIES)


def test_genetic_code_integrity():
    assert len(CODON_TABLE) == 64
    assert CODON_TABLE["ATG"] == "M"
    assert CODON_TABLE["TGG"] == "W"
    assert CODON_TABLE["TAA"] == "*"
    assert CODON_TABLE["CTG"] == "L"
    assert CODON_TABLE["AAA"] == "K"
    assert CODON_TABLE["GCT"] == "A"
    # degeneracy classes of the standard code: 2x1, 9x2, 1x3, 5x4, 3x6
    sizes = {k: len(v) for k, v in DEGENERACY_CLASSES.items()}
    assert sizes == {1: 2, 2: 9, 3: 1, 4: 5, 6: 3}
    assert sum(len(v) for v in FAMILIES.values()) == 61


def test_counting_rules():
    out = analysis.count_codons("ATGGCTGCTGCTTAA")
    assert out["counts"] == {"ATG": 1, "GCT": 3, "TAA": 1}
    assert out["stop_codons"] == 1 and out["leftover_bases"] == 0
    out = analysis.count_codons("ATGNNNGCT")
    assert out["skipped_codons"] == 1
    assert out["counts"] == {"ATG": 1, "GCT": 1}
    out = analysis.count_codons("ATGGC")
    assert out["leftover_bases"] == 2 and out["counts"] == {"ATG": 1}


def test_rscu_hand_computed():
    # Ala family: GCT x3, nothing else -> RSCU(GCT) = 4 * 3/3 = 4.0
    counts = {"ATG": 1, "GCT": 3}
    r = analysis.rscu(counts)
    assert r["GCT"] == 4.0
    assert r["GCC"] == 0.0 and r["GCA"] == 0.0 and r["GCG"] == 0.0
    assert r["ATG"] == 1.0  # single-codon family
    # equal use of a 2-fold family -> both RSCU exactly 1
    r2 = analysis.rscu({"AAA": 7, "AAG": 7})
    assert r2["AAA"] == 1.0 and r2["AAG"] == 1.0
    # fractions and per-thousand agree with the same counts
    fr = analysis.family_fraction(counts)
    assert fr["GCT"] == 1.0
    pt = analysis.per_thousand(counts)
    assert pt["GCT"] == 750.0 and pt["ATG"] == 250.0


def test_gc_positions_hand_computed():
    # codons ATG GCT GCC: GC bases 1+2+3=6 of 9; thirds G,T,C -> 2/3
    bundle = analysis.analyse_sequence("ATGGCTGCC")
    gc = bundle["gc"]
    assert abs(gc["gc"] - 6 / 9) < 1e-12
    assert abs(gc["gc1"] - 2 / 3) < 1e-12
    assert abs(gc["gc2"] - 2 / 3) < 1e-12
    assert abs(gc["gc3"] - 2 / 3) < 1e-12


def test_nc_extremes_and_fallback():
    # one codon per amino acid only -> maximal bias -> Nc = 20
    one_each = {}
    for aa, codons in FAMILIES.items():
        one_each[codons[0]] = 10
    assert abs(analysis.effective_number_of_codons(one_each) - 20.0) < 1e-9
    # perfectly uniform usage -> Nc hits the 61 ceiling
    uniform = {c: 100 for codons in FAMILIES.values() for c in codons}
    assert analysis.effective_number_of_codons(uniform) == 61.0
    # only Lys (AAA/AAG 5/5): F = (10*0.5 - 1)/9 = 4/9; missing classes
    # fall back to that mean -> Nc = 2 + 18/(4/9) = 42.5
    nc = analysis.effective_number_of_codons({"AAA": 5, "AAG": 5})
    assert abs(nc - 42.5) < 1e-9
    # nothing usable at all
    assert analysis.effective_number_of_codons({}) is None
    assert analysis.effective_number_of_codons({"ATG": 1}) is None


if __name__ == "__main__":
    test_genetic_code_integrity()
    test_counting_rules()
    test_rscu_hand_computed()
    test_gc_positions_hand_computed()
    test_nc_extremes_and_fallback()
    print("test_analysis: all passed")
