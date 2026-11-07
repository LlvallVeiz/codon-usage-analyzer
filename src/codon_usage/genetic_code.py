"""The standard genetic code (NCBI translation table 1).

The 64 codons are generated in the canonical order (bases T, C, A, G;
first position slowest) and zipped with the equally canonical amino
acid string, exactly as published by NCBI [see docs/RESEARCH.md].
Stop codons map to "*". `CODON_TABLE` maps codon -> one-letter amino
acid; `FAMILIES` groups sense codons by amino acid, which is the
structure every codon-usage statistic in this project is built on.
"""

from __future__ import annotations

_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*W"   # T..
        "LLLLPPPPHHQQRRRR"   # C..
        "IIIMTTTTNNKKSSRR"   # A..
        "VVVVAAAADDEEGGGG")  # G..

CODONS: list[str] = [a + b + c for a in _BASES for b in _BASES
                     for c in _BASES]
CODON_TABLE: dict[str, str] = dict(zip(CODONS, _AAS))

STOP_CODONS: list[str] = [c for c in CODONS if CODON_TABLE[c] == "*"]

# Amino acid (one-letter) -> its synonymous sense codons.
FAMILIES: dict[str, list[str]] = {}
for _codon in CODONS:
    _aa = CODON_TABLE[_codon]
    if _aa != "*":
        FAMILIES.setdefault(_aa, []).append(_codon)

# Degeneracy classes used by Wright's Nc: {degeneracy: [amino acids]}.
DEGENERACY_CLASSES: dict[int, list[str]] = {}
for _aa, _codons in FAMILIES.items():
    DEGENERACY_CLASSES.setdefault(len(_codons), []).append(_aa)

AMINO_ACID_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}
