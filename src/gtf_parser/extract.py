"""Extract custom fields from parsed GTF/GFF3 records, with deduplication."""

from __future__ import annotations

from collections.abc import Mapping, Set
import sys
from pathlib import Path
from typing import TextIO

from .parser import parse


def extract(
    source: str | Path | TextIO,
    fields: list[str],
    feature_types: set[str] | None = None,
    output: TextIO | None = None,
    separator: str = "\t",
    deduplicate: bool = True,
    index_path: str | Path | None = None,
    attribute_filters: Mapping[str, Set[str]] | None = None,
) -> None:
    """Extract selected fields from a GTF/GFF3 file and write them out.

    Parameters
    ----------
    source:
        Input GTF/GFF3 file path or open stream.
    fields:
        List of column names or attribute keys to extract.
    feature_types:
        Restrict to these feature types (``None`` = all).
    output:
        Writable text stream for output (default: stdout).
    separator:
        Column separator for the output.
    deduplicate:
        If ``True``, skip duplicate rows.
    index_path:
        Optional SQLite index path created with ``gtf-parser index``.
    attribute_filters:
        Optional attribute filters in the form ``{key: {value1, value2}}``.
    """
    out = output or sys.stdout
    seen: set[tuple[str, ...]] | None = set() if deduplicate else None

    out.write(separator.join(fields) + "\n")

    for record in parse(
        source,
        feature_types=feature_types,
        index_path=index_path,
        attribute_filters=attribute_filters,
    ):
        values = tuple(record.get(f) or "" for f in fields)
        if seen is not None:
            if values in seen:
                continue
            seen.add(values)
        out.write(separator.join(values) + "\n")
