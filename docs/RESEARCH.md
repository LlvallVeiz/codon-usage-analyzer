# Background & Research — Codon Usage Analyzer

This document records the primary sources every statistic in this
project is built on, the exact formulas implemented, and the places
where the code makes a documented choice between conventions in the
literature. Every source below was consulted while writing the code.

## 1. The genetic code: NCBI translation table 1

All codon → amino-acid assignments use the **standard genetic code**,
NCBI translation table 1 (`transl_table=1`), as published in the NCBI
Taxonomy "Genetic Codes" reference [1]. The implementation stores the
table in the same canonical compact form NCBI publishes it in: codons
ordered with bases T, C, A, G (first position slowest) zipped against
the amino-acid string

```text
FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
```

A unit test re-derives the degeneracy structure from the table and
asserts the textbook values — 2 one-fold (Met, Trp), 9 two-fold,
1 three-fold (Ile), 5 four-fold and 3 six-fold amino acids, 61 sense
codons — so a typo in the table cannot slip through silently.

## 2. RSCU: Sharp & Li (1987), after Sharp, Tuohy & Mosurski (1986)

The **relative synonymous codon usage** of a codon is "the observed
frequency of the codon divided by the frequency expected under the
assumption of equal usage of the synonymous codons for an amino
acid" — the definition as stated by Sharp & Li [2], who took the
statistic from Sharp, Tuohy & Mosurski's yeast analysis [3]:

```text
RSCU_ij  =  X_ij / ( (1/n_i) · Σ_j X_ij )  =  n_i · X_ij / Σ_j X_ij
```

where X_ij is the count of the j-th codon for amino acid i and n_i is
that amino acid's degeneracy. Properties the implementation (and its
tests) rely on: RSCU = 1 for every codon under perfectly uniform
synonymous usage; RSCU sums to n_i within each family; the maximum
for one codon is n_i (it takes the whole family). RSCU is computed
here for the 61 sense codons; stop codons are counted but excluded,
as is conventional — they are not translated, and pooling them into a
"family" would mix termination biology into a translation statistic.
Sharp & Li's paper is equally the source of the Codon Adaptation
Index (CAI); CAI is deliberately *not* implemented, because it needs
a reference set of highly expressed genes — shipping it without one
would produce a number with no defensible meaning (see README,
"Limitations").

## 3. Nc: Wright (1990)

The **effective number of codons** "quantifies how far the codon
usage of a gene departs from equal usage of synonymous codons… Nc can
take values from 20, in the case of extreme bias where one codon is
exclusively used for each aa, to 61 when the use of alternative
synonymous codons is equally likely" [4]. Wright's estimator:

```text
Nc  =  2 + 9/F̄2 + 1/F̄3 + 5/F̄4 + 3/F̄6
F   =  ( n · Σ p_i² − 1 ) / ( n − 1 )        (per amino-acid family)
```

where F̄k is the mean homozygosity F over the amino acids with
k-fold degeneracy (9, 1, 5 and 3 amino acids respectively; the
leading 2 accounts for Met and Trp), p_i are the within-family codon
proportions, and n is the family's total count. The (n−1) correction
is Wright's sampling-without-replacement estimator; a later
resampling study confirmed it behaves better at short gene lengths
than the with-replacement variant [5], and improved estimators have
been proposed since [6] — this project implements Wright's original,
clamped to [20, 61].

Documented choices where Wright's paper leaves implementation room:

* **Families with n < 2** cannot produce F (division by n−1) and are
  excluded from their class mean.
* **A degeneracy class with no usable family** takes the mean of the
  available class means as its F̄. The alternative (dropping the
  term) would silently change the statistic's definition per gene;
  the fallback keeps Nc comparable across genes of different
  composition, at the cost of some accuracy for very short inputs —
  Nc on a handful of codons is a rough indicator, and the README
  says so.
* **F is not floored at 1/k.** The corrected estimator can fall
  below the uniform-usage value through sampling noise; flooring
  would bias Nc upward for short genes. The one true singularity,
  a class mean of exactly 0 (every observed codon seen at most
  once), is handled by letting that class contribute its maximum,
  n_families × degeneracy, before the [20, 61] clamp.

The committed tests pin all three behaviours to hand-computed
values: single-codon-per-family counts give Nc = 20 exactly; uniform
counts give the 61 ceiling; a lone Lys family split 5/5 gives
F = 4/9 and Nc = 42.5 through the fallback path.

## 4. GC content and GC3

Overall GC and the positional GC1/GC2/GC3 are simple (G+C) fractions;
GC3 — the GC fraction at third codon positions — is the quantity
Wright's Nc-plot is drawn against [4], because third positions are
where synonymous choice lives. Here GC3 is computed over the third
positions of **all sense codons** (stop codons excluded). Some tools
restrict GC3 to synonymous ("GC3s") sites only; for sequences whose
Met/Trp content is small the two are close, but they are not the
same statistic, and this project's reports label plainly which one
is printed.

## 5. The sample data is synthetic — and labelled as such

`examples/sample.fasta` contains two artificial CDSs built for this
project: sequence A encodes a 20-aa block three times using one
strongly preferred codon per amino acid (Nc = 20.00, the bias
extreme), sequence B encodes a different block with the *opposite*
preferences. Pooled, the biases partially cancel (Nc = 35.90) — a
compact, checkable demonstration that Nc responds to usage, not to
amino-acid composition. Every number was cross-checked against an
independent throwaway reimplementation before being committed.

## References

1. NCBI Taxonomy, *The Genetic Codes* (translation table 1) —
   https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
2. Sharp PM, Li WH (1987). The codon adaptation index — a measure of
   directional synonymous codon usage bias, and its potential
   applications. *Nucleic Acids Research* 15(3):1281–1295.
   DOI 10.1093/nar/15.3.1281 —
   https://pubmed.ncbi.nlm.nih.gov/3547335/
3. Sharp PM, Tuohy TMF, Mosurski KR (1986). Codon usage in yeast:
   cluster analysis clearly differentiates highly and lowly expressed
   genes. *Nucleic Acids Research* 14(13):5125–5143.
   DOI 10.1093/nar/14.13.5125
4. Wright F (1990). The 'effective number of codons' used in a gene.
   *Gene* 87(1):23–29. DOI 10.1016/0378-1119(90)90491-9 —
   https://pubmed.ncbi.nlm.nih.gov/2110097/
5. Fuglsang A (2006). Estimating the "effective number of codons":
   the Wright way of determining codon homozygosity leads to superior
   estimates. *Genetics* 172(2):1301–1307.
   DOI 10.1534/genetics.105.049643 —
   https://pubmed.ncbi.nlm.nih.gov/16299393/
6. Sun X, Yang Q, Xia X (2013). An improved implementation of the
   effective number of codons (Nc). *Molecular Biology and Evolution*
   30(1):191–196. DOI 10.1093/molbev/mss201
