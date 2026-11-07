"""A small, strict-enough FASTA parser (standard library only).

Records are returned as (header, sequence) with the sequence
upper-cased and whitespace removed. The header is the full line after
">"; the record id used in reports is its first whitespace-delimited
token, matching common practice.
"""

from __future__ import annotations


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into [(header, sequence), ...]."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:]
            chunks = []
        else:
            if header is None:
                raise ValueError("sequence data before first '>' header")
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def read_fasta(path: str) -> list[tuple[str, str]]:
    with open(path) as fh:
        return parse_fasta(fh.read())


def record_id(header: str) -> str:
    return header.split()[0] if header.split() else header
