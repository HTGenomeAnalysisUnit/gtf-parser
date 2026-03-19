"""Core GTF/GFF3 parser.

Handles both GTF (Gene Transfer Format) and GFF3 (General Feature Format v3)
files, auto-detecting the format based on attribute syntax.

GTF attributes:  key "value"; key2 "value2";
GFF3 attributes: key=value;key2=value2
"""

from __future__ import annotations

import gzip
import json
import os
import re
import sqlite3
from collections.abc import Mapping, Set
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Iterator, TextIO

# The 9 standard columns in GTF/GFF files
COLUMNS = ("seqname", "source", "feature", "start", "end", "score", "strand", "frame")
INDEX_SCHEMA_VERSION = "2"
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
_GTF_DETECT_RE = re.compile(r'"\s*;')


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
    for match in _GTF_ATTR_RE.finditer(raw):
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
    if "=" in attr_string and not _GTF_DETECT_RE.search(attr_string):
        return "gff3"
    return "gtf"


def _open_file(path: str | Path) -> IO[str]:
    """Open a plain or gzip-compressed file."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def default_index_path(source: str | Path) -> Path:
    """Return the default on-disk index path for an input annotation file."""
    source = Path(source)
    return Path(f"{source}.idx.sqlite")


def _source_fingerprint(path: Path) -> tuple[str, str]:
    stat = path.stat()
    return str(stat.st_size), str(stat.st_mtime_ns)


def _parse_meta_set(raw: str | None) -> set[str] | None:
    """Parse comma-separated metadata sets; ``None`` means unconstrained."""
    if raw is None or raw == "*":
        return None
    return {item for item in raw.split(",") if item}


def _index_status(
    source: Path,
    index_path: Path,
    feature_types: set[str] | None = None,
    attribute_filters: dict[str, set[str]] | None = None,
) -> tuple[bool, str]:
    if not index_path.exists():
        return False, "index file does not exist"
    if not source.exists():
        return False, "source file does not exist"

    try:
        conn = sqlite3.connect(index_path)
        try:
            rows = conn.execute("SELECT key, value FROM meta").fetchall()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        return False, f"cannot read index metadata: {exc}"

    meta = {key: value for key, value in rows}
    schema = meta.get("schema_version")
    if schema != INDEX_SCHEMA_VERSION:
        return False, f"schema mismatch (index={schema}, expected={INDEX_SCHEMA_VERSION})"

    source_size, source_mtime_ns = _source_fingerprint(source)
    if meta.get("source_path") != str(source.resolve()):
        return False, "index was built for a different source path"
    if meta.get("source_size") != source_size or meta.get("source_mtime_ns") != source_mtime_ns:
        return False, "index is stale (source file changed since index build)"

    indexed_features = _parse_meta_set(meta.get("indexed_feature_types"))
    if feature_types and indexed_features is not None:
        missing_features = sorted(feature_types - indexed_features)
        if missing_features:
            return (
                False,
                "index missing feature type(s): " + ", ".join(missing_features),
            )

    indexed_attr_keys = _parse_meta_set(meta.get("indexed_attribute_keys"))
    if attribute_filters and indexed_attr_keys is not None:
        missing_keys = sorted(set(attribute_filters) - indexed_attr_keys)
        if missing_keys:
            return (
                False,
                "index missing attribute key(s): " + ", ".join(missing_keys),
            )

    return True, "ok"


def index_status(
    source: str | Path,
    index_path: str | Path,
    feature_types: set[str] | None = None,
    attribute_filters: Mapping[str, Set[str]] | None = None,
) -> tuple[bool, str]:
    """Validate whether *index_path* can be used for *source*."""
    normalized_filters = _normalize_attribute_filters(attribute_filters)
    return _index_status(
        Path(source),
        Path(index_path),
        feature_types=feature_types,
        attribute_filters=normalized_filters,
    )


def _normalize_attribute_filters(
    attribute_filters: Mapping[str, Set[str]] | None,
) -> dict[str, set[str]] | None:
    if not attribute_filters:
        return None

    normalized: dict[str, set[str]] = {}
    for key, values in attribute_filters.items():
        clean_key = key.strip()
        if not clean_key:
            continue

        clean_values = {value.strip() for value in values if value.strip()}
        if not clean_values:
            continue

        if clean_key in normalized:
            normalized[clean_key].update(clean_values)
        else:
            normalized[clean_key] = clean_values

    return normalized or None


def _record_matches_attribute_filters(
    attributes: dict[str, str],
    attribute_filters: dict[str, set[str]] | None,
) -> bool:
    if not attribute_filters:
        return True
    for key, allowed_values in attribute_filters.items():
        if attributes.get(key) not in allowed_values:
            return False
    return True


def _parse_stream(
    fh: TextIO,
    feature_types: set[str] | None = None,
    attribute_filters: dict[str, set[str]] | None = None,
) -> Iterator[Record]:
    """Parse an already-open GTF/GFF3 text stream."""
    fmt: str | None = None  # auto-detect on first data line
    parse_attrs = _parse_gtf_attributes

    for line in fh:
        line = line.rstrip("\n\r")
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t", 8)
        if len(parts) < 9:
            continue

        feature = parts[2]
        if feature_types and feature not in feature_types:
            continue

        attr_raw = parts[8]

        if fmt is None:
            fmt = _detect_format(attr_raw)
            parse_attrs = _parse_gff3_attributes if fmt == "gff3" else _parse_gtf_attributes

        attributes = parse_attrs(attr_raw)
        if not _record_matches_attribute_filters(attributes, attribute_filters):
            continue

        yield Record(
            seqname=parts[0],
            source=parts[1],
            feature=feature,
            start=int(parts[3]),
            end=int(parts[4]),
            score=parts[5],
            strand=parts[6],
            frame=parts[7],
            attributes=attributes,
        )


def _index_matches_source(
    source: Path,
    index_path: Path,
    feature_types: set[str] | None = None,
    attribute_filters: dict[str, set[str]] | None = None,
) -> bool:
    usable, _ = _index_status(
        source,
        index_path,
        feature_types=feature_types,
        attribute_filters=attribute_filters,
    )
    return usable


def _iter_from_index(
    index_path: Path,
    feature_types: set[str] | None = None,
    attribute_filters: dict[str, set[str]] | None = None,
) -> Iterator[Record]:
    conn = sqlite3.connect(index_path)
    try:
        clauses: list[str] = []
        params: list[str] = []

        if feature_types:
            ordered_features = sorted(feature_types)
            placeholders = ",".join("?" for _ in ordered_features)
            clauses.append(f"r.feature IN ({placeholders})")
            params.extend(ordered_features)

        if attribute_filters:
            subqueries: list[str] = []
            subquery_params: list[str] = []
            for key, values in sorted(attribute_filters.items()):
                ordered_values = sorted(values)
                placeholders = ",".join("?" for _ in ordered_values)
                subqueries.append(
                    "SELECT row_num FROM attributes "
                    f"WHERE attr_key = ? AND attr_value IN ({placeholders})"
                )
                subquery_params.append(key)
                subquery_params.extend(ordered_values)

            intersected = " INTERSECT ".join(subqueries)
            clauses.append(f"r.row_num IN ({intersected})")
            params.extend(subquery_params)

        where_clause = f" WHERE {' AND '.join(clauses)}" if clauses else ""
        query = (
            "SELECT r.seqname, r.source, r.feature, r.start, r.end, r.score, r.strand, r.frame, r.attributes_json "
            "FROM records r"
            f"{where_clause} ORDER BY r.row_num"
        )
        cursor = conn.execute(query, params)

        for row in cursor:
            yield Record(
                seqname=row[0],
                source=row[1],
                feature=row[2],
                start=row[3],
                end=row[4],
                score=row[5],
                strand=row[6],
                frame=row[7],
                attributes=json.loads(row[8]),
            )
    finally:
        conn.close()


def build_index(
    source: str | Path,
    index_path: str | Path | None = None,
    force: bool = False,
    feature_types: set[str] | None = None,
    attribute_keys: set[str] | None = None,
) -> Path:
    """Build an SQLite index for fast retrieval.

    Parameters
    ----------
    source:
        Input GTF/GFF3 file path.
    index_path:
        Output SQLite path. Defaults to ``<source>.idx.sqlite``.
    force:
        Overwrite the output index if it already exists.
    feature_types:
        Optional subset of feature types to index.
    attribute_keys:
        Optional subset of attribute keys to index for ``attribute_filters``.
    """
    source_path = Path(source)
    if not source_path.exists():
        raise FileNotFoundError(source_path)

    index = Path(index_path) if index_path is not None else default_index_path(source_path)
    if index.exists():
        if force:
            index.unlink()
        else:
            raise FileExistsError(index)

    index.parent.mkdir(parents=True, exist_ok=True)
    tmp_index = index.with_suffix(index.suffix + ".tmp")
    if tmp_index.exists():
        tmp_index.unlink()

    source_size, source_mtime_ns = _source_fingerprint(source_path)
    normalized_attribute_keys = (
        {key.strip() for key in attribute_keys if key.strip()} if attribute_keys else None
    )

    conn: sqlite3.Connection | None = None
    try:
        conn = sqlite3.connect(tmp_index)
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA temp_store=MEMORY")
        conn.execute(
            "CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        conn.execute(
            "CREATE TABLE records ("
            "row_num INTEGER PRIMARY KEY, "
            "seqname TEXT NOT NULL, "
            "source TEXT NOT NULL, "
            "feature TEXT NOT NULL, "
            "start INTEGER NOT NULL, "
            "end INTEGER NOT NULL, "
            "score TEXT NOT NULL, "
            "strand TEXT NOT NULL, "
            "frame TEXT NOT NULL, "
            "attributes_json TEXT NOT NULL)"
        )
        conn.execute(
            "CREATE TABLE attributes ("
            "row_num INTEGER NOT NULL, "
            "attr_key TEXT NOT NULL, "
            "attr_value TEXT NOT NULL)"
        )
        conn.execute("CREATE INDEX idx_records_feature ON records(feature)")
        conn.execute(
            "CREATE INDEX idx_attributes_key_value_row "
            "ON attributes(attr_key, attr_value, row_num)"
        )

        conn.executemany(
            "INSERT INTO meta(key, value) VALUES(?, ?)",
            (
                ("schema_version", INDEX_SCHEMA_VERSION),
                ("source_path", str(source_path.resolve())),
                ("source_size", source_size),
                ("source_mtime_ns", source_mtime_ns),
                (
                    "indexed_feature_types",
                    "*" if feature_types is None else ",".join(sorted(feature_types)),
                ),
                (
                    "indexed_attribute_keys",
                    "*"
                    if normalized_attribute_keys is None
                    else ",".join(sorted(normalized_attribute_keys)),
                ),
            ),
        )

        row_num = 0
        batch: list[tuple[int, str, str, str, int, int, str, str, str, str]] = []
        attr_batch: list[tuple[int, str, str]] = []
        with _open_file(source_path) as fh:
            for record in _parse_stream(fh, feature_types=feature_types):
                row_num += 1
                batch.append(
                    (
                        row_num,
                        record.seqname,
                        record.source,
                        record.feature,
                        record.start,
                        record.end,
                        record.score,
                        record.strand,
                        record.frame,
                        json.dumps(record.attributes, separators=(",", ":")),
                    )
                )
                if normalized_attribute_keys is None:
                    for key, value in record.attributes.items():
                        attr_batch.append((row_num, key, value))
                else:
                    for key in normalized_attribute_keys:
                        value = record.attributes.get(key)
                        if value is not None:
                            attr_batch.append((row_num, key, value))

                if len(batch) >= 10_000:
                    conn.executemany(
                        "INSERT INTO records("
                        "row_num, seqname, source, feature, start, end, score, strand, frame, attributes_json"
                        ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    batch.clear()
                if len(attr_batch) >= 50_000:
                    conn.executemany(
                        "INSERT INTO attributes(row_num, attr_key, attr_value) VALUES (?, ?, ?)",
                        attr_batch,
                    )
                    attr_batch.clear()

        if batch:
            conn.executemany(
                "INSERT INTO records("
                "row_num, seqname, source, feature, start, end, score, strand, frame, attributes_json"
                ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                batch,
            )
        if attr_batch:
            conn.executemany(
                "INSERT INTO attributes(row_num, attr_key, attr_value) VALUES (?, ?, ?)",
                attr_batch,
            )

        conn.commit()
    except Exception:
        if conn is not None:
            conn.close()
        if tmp_index.exists():
            tmp_index.unlink()
        raise
    else:
        conn.close()

    os.replace(tmp_index, index)
    return index


def parse(
    source: str | Path | TextIO,
    feature_types: set[str] | None = None,
    index_path: str | Path | None = None,
    attribute_filters: Mapping[str, Set[str]] | None = None,
) -> Iterator[Record]:
    """Parse a GTF or GFF3 file, yielding Record objects.

    Parameters
    ----------
    source:
        File path (str/Path) or an already-open text stream.
    feature_types:
        If provided, only records whose ``feature`` column matches one of the
        given types are yielded.  Pass ``None`` to yield all records.
    index_path:
        Optional path to an SQLite index created by :func:`build_index`.
        If provided and valid for *source*, parsing is served from the index.
    attribute_filters:
        Optional attribute filters in the form ``{key: {value1, value2}}``.
        A record matches only if all keys match one of their allowed values.
    """
    normalized_filters = _normalize_attribute_filters(attribute_filters)

    if isinstance(source, (str, Path)):
        source_path = Path(source)
        if index_path is not None:
            index = Path(index_path)
            if _index_matches_source(
                source_path,
                index,
                feature_types=feature_types,
                attribute_filters=normalized_filters,
            ):
                yield from _iter_from_index(
                    index,
                    feature_types=feature_types,
                    attribute_filters=normalized_filters,
                )
                return

        fh = _open_file(source_path)
        should_close = True
    else:
        fh = source
        should_close = False

    try:
        yield from _parse_stream(
            fh,
            feature_types=feature_types,
            attribute_filters=normalized_filters,
        )
    finally:
        if should_close:
            fh.close()
