"""Convert GTF/GFF3 records to BED format with configurable ID column."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TextIO

from .parser import Record, parse


def to_bed(
    source: str | Path | TextIO,
    feature_type: str,
    id_fields: list[str],
    id_separator: str = "|",
    output: TextIO | None = None,
    deduplicate: bool = True,
) -> None:
    """Convert GTF/GFF3 records of a given feature type to BED6 format.

    BED columns produced:
        chrom, chromStart (0-based), chromEnd, name, score, strand

    The ``name`` column is built by joining the values of *id_fields* (which
    can reference main columns or attribute keys) with *id_separator*.

    Parameters
    ----------
    source:
        Input file path or open stream.
    feature_type:
        Only records of this feature type are converted.
    id_fields:
        Attribute keys / column names whose values are concatenated to form
        the BED ``name`` column.
    id_separator:
        String used to join id_fields values.
    output:
        Writable text stream (default: stdout).
    deduplicate:
        If ``True``, skip duplicate BED lines.
    """
    out = output or sys.stdout
    seen: set[str] | None = set() if deduplicate else None

    for record in parse(source, feature_types={feature_type}):
        name = id_separator.join(record.get(f) or "" for f in id_fields)
        # BED is 0-based half-open; GTF/GFF are 1-based inclusive
        bed_start = record.start - 1
        bed_end = record.end
        score = record.score if record.score != "." else "0"
        line = f"{record.seqname}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{record.strand}"

        if seen is not None:
            if line in seen:
                continue
            seen.add(line)

        out.write(line + "\n")
