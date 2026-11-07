# Codon Usage Analyzer

Codon-usage analysis for coding sequences in **pure Python 3
(standard library only)**: per-sequence and pooled codon counts, a
full codon-usage table (count, per-thousand, within-family fraction,
**RSCU**), **GC / GC1 / GC2 / GC3** content, and **Wright's Nc**
(effective number of codons) — with a text bar chart of RSCU values.

## Background

Genes do not use synonymous codons equally. Two classical statistics
quantify the bias (full derivations and sources in
[docs/RESEARCH.md](docs/RESEARCH.md)):

**RSCU** — relative synonymous codon usage (Sharp & Li 1987):

```text
RSCU_ij = n_i · X_ij / Σ_j X_ij
```

the observed count of a codon divided by the count expected if all
synonyms of its amino acid were used equally. RSCU = 1 means no
preference; the maximum is the family's degeneracy (one codon takes
the whole family).

**Nc** — the effective number of codons (Wright 1990):

```text
Nc = 2 + 9/F̄2 + 1/F̄3 + 5/F̄4 + 3/F̄6,
F  = ( n·Σp² − 1 ) / ( n − 1 )   per amino-acid family
```

Nc runs from **20** (extreme bias: one codon per amino acid) to
**61** (no bias at all), independent of gene length and amino-acid
composition — which is what makes it comparable across genes.

Counting conventions, stated plainly: sequences are read in frame 0
as CDS; complete codons only (trailing partial bases are reported);
codons with ambiguous bases are skipped and reported; stop codons
are counted but excluded from RSCU, GC3 and Nc. Codon assignments
are the standard genetic code (NCBI translation table 1).

## Project structure

```text
codon-usage-analyzer/
├── analyze.py              # entry point
├── src/codon_usage/
│   ├── genetic_code.py     # NCBI table 1 + synonymous families
│   ├── fasta.py            # FASTA parser
│   ├── analysis.py         # counts, RSCU, GC/GC3, Wright's Nc
│   ├── report.py           # table (text/CSV) + RSCU bar chart
│   └── cli.py              # argparse wiring
├── tests/
│   └── test_analysis.py    # hand-computed expectations throughout
├── examples/
│   ├── sample.fasta        # two synthetic CDSs, opposite biases
│   ├── codon_table.txt     # generated reports (pooled)
│   ├── codon_table.csv
│   ├── rscu_chart.txt
│   └── example_output.txt  # captured run
└── docs/RESEARCH.md        # sources, formulas, documented choices
```

## How to run

```bash
python3 analyze.py examples/sample.fasta --outdir examples
python3 tests/test_analysis.py
```

## Example output

The sample file holds two *synthetic* CDSs, each built to use one
preferred codon per amino acid — sequence A preferring one set,
sequence B the opposite set:

```text
--- synthetic_bias_A ---
sense codons : 60     GC3 : 0.250     Nc (Wright) : 20.00
--- synthetic_bias_B ---
sense codons : 20     GC3 : 0.850     Nc (Wright) : 20.00
--- GLOBAL (pooled, 2 sequences) ---
sense codons : 80     GC3 : 0.400     Nc (Wright) : 35.90
```

Each sequence alone sits at the Nc floor (20 = one codon per amino
acid, by construction), while the pooled analysis lands mid-scale —
the two opposite biases partially cancel. It is a compact proof that
the statistic responds to *usage*, not composition: in sequence A,
leucine is 6/6 CTG (RSCU 6.00, the maximum), and every one of these
numbers was cross-checked against an independent reimplementation
before being committed.

## Design decisions

- **Hand-computed tests, not snapshots** — the suite asserts values
  worked out by hand (RSCU 4.00 for a 3-of-3 Ala family, GC3 = 2/3
  for ATG/GCT/GCC, Nc = 20 / 61 / 42.5 for constructed count
  profiles), so a formula error fails loudly instead of being
  ratified by its own output.
- **Documented Nc edge behaviour** — families with < 2 observations
  are excluded from their class mean; an entirely missing class
  falls back to the mean of the available classes; F is not floored
  at 1/k (small-sample honesty over cosmetic stability). Rationale
  and trade-offs: docs/RESEARCH.md §3.
- **No CAI** — Sharp & Li's Codon Adaptation Index needs a reference
  set of highly expressed genes; without one the number would be
  meaningless, so it is left out on purpose rather than faked.

## Limitations & future work

- Input is assumed to be in-frame CDS; the tool does not search for
  ORFs or validate start/stop placement (Day-01's
  `dna-sequence-toolkit` covers ORF finding).
- Only the standard genetic code is supported; mitochondrial and
  other NCBI tables would slot into `genetic_code.py` cleanly.
- Nc for very short sequences (< ~50 codons) is inherently noisy —
  Wright's estimator is unbiased-ish, not magic.
- Natural extensions: per-sequence CSV export, an Nc-vs-GC3 plot
  (Wright's Nc-plot), and CAI with a user-supplied reference set.
