"""Core GTF/GFF3 parser.

Handles both GTF (Gene Transfer Format) and GFF3 (General Feature Format v3)
files, auto-detecting the format based on attribute syntax.

GTF attributes:  key "value"; key2 "value2";
GFF3 attributes: key=value;key2=value2
"""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Iterator, TextIO

# The 9 standard columns in GTF/GFF files
COLUMNS = ("seqname", "source", "feature", "start", "end", "score", "strand", "frame")


@dataclass(slots=True)
class Record:
    """A single GTF/GFF3 record (one line)."""

    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attributes: dict[str, str] = field(default_factory=dict)

    def get(self, key: str) -> str | None:
        """Get a value from main columns or attributes."""
        if key in COLUMNS:
            return str(getattr(self, key))
        return self.attributes.get(key)


def _parse_gtf_attributes(raw: str) -> dict[str, str]:
    """Parse GTF-style attributes: key "value"; key2 "value2";"""
    attrs: dict[str, str] = {}
    for match in re.finditer(r'(\w+)\s+"([^"]*)"', raw):
        attrs[match.group(1)] = match.group(2)
    return attrs


def _parse_gff3_attributes(raw: str) -> dict[str, str]:
    """Parse GFF3-style attributes: key=value;key2=value2"""
    attrs: dict[str, str] = {}
    for part in raw.split(";"):
        part = part.strip()
        if "=" in part:
            key, _, value = part.partition("=")
            attrs[key.strip()] = value.strip()
    return attrs


def _detect_format(attr_string: str) -> str:
    """Detect whether an attribute string is GTF or GFF3 format."""
    if "=" in attr_string and not re.search(r'"\s*;', attr_string):
        return "gff3"
    return "gtf"


def _open_file(path: str | Path) -> IO[str]:
    """Open a plain or gzip-compressed file."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def parse(
    source: str | Path | TextIO,
    feature_types: set[str] | None = None,
) -> Iterator[Record]:
    """Parse a GTF or GFF3 file, yielding Record objects.

    Parameters
    ----------
    source:
        File path (str/Path) or an already-open text stream.
    feature_types:
        If provided, only records whose ``feature`` column matches one of the
        given types are yielded.  Pass ``None`` to yield all records.
    """
    if isinstance(source, (str, Path)):
        fh = _open_file(source)
        should_close = True
    else:
        fh = source
        should_close = False

    fmt: str | None = None  # auto-detect on first data line

    try:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            feature = parts[2]
            if feature_types and feature not in feature_types:
                continue

            attr_raw = parts[8]

            if fmt is None:
                fmt = _detect_format(attr_raw)
            parse_attrs = (
                _parse_gff3_attributes if fmt == "gff3" else _parse_gtf_attributes
            )

            yield Record(
                seqname=parts[0],
                source=parts[1],
                feature=feature,
                start=int(parts[3]),
                end=int(parts[4]),
                score=parts[5],
                strand=parts[6],
                frame=parts[7],
                attributes=parse_attrs(attr_raw),
            )
    finally:
        if should_close:
            fh.close()
