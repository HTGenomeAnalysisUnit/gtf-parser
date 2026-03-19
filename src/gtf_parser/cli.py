"""Command-line interface for gtf-parser."""

from __future__ import annotations

import argparse
import sys

from . import __version__


def _add_common_args(sub: argparse.ArgumentParser) -> None:
    """Arguments shared across sub-commands."""
    sub.add_argument("input", help="Input GTF or GFF3 file (plain or .gz)")
    sub.add_argument(
        "-t",
        "--feature-types",
        nargs="+",
        default=None,
        help="Feature type(s) to process (default: all)",
    )
    sub.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output file (default: stdout)",
    )
    sub.add_argument(
        "--no-dedup",
        action="store_true",
        help="Disable deduplication of output rows",
    )
    sub.add_argument(
        "--index",
        default=None,
        help="Optional SQLite index path created with `gtf-parser index`",
    )
    sub.add_argument(
        "--where",
        action="append",
        default=None,
        help=(
            "Attribute filter in the form key=value or key=v1,v2. "
            "Can be repeated to combine filters."
        ),
    )


def _parse_where_filters(where_filters: list[str] | None) -> dict[str, set[str]] | None:
    if not where_filters:
        return None

    parsed: dict[str, set[str]] = {}
    for raw_filter in where_filters:
        key, sep, raw_values = raw_filter.partition("=")
        if not sep:
            raise ValueError(f"Invalid --where value '{raw_filter}': expected key=value")

        clean_key = key.strip()
        if not clean_key:
            raise ValueError(f"Invalid --where value '{raw_filter}': empty key")

        values = {value.strip() for value in raw_values.split(",") if value.strip()}
        if not values:
            raise ValueError(f"Invalid --where value '{raw_filter}': empty value list")

        parsed.setdefault(clean_key, set()).update(values)

    return parsed


def _warn_if_index_unusable(
    source: str,
    index_path: str | None,
    feature_types: set[str] | None = None,
    attribute_filters: dict[str, set[str]] | None = None,
) -> None:
    if not index_path:
        return

    from .parser import index_status

    usable, reason = index_status(
        source,
        index_path,
        feature_types=feature_types,
        attribute_filters=attribute_filters,
    )
    if not usable:
        print(
            f"Warning: index '{index_path}' will be ignored ({reason}).",
            file=sys.stderr,
        )


def _run_extract(args: argparse.Namespace) -> None:
    from .extract import extract

    try:
        where_filters = _parse_where_filters(args.where)
    except ValueError as exc:
        raise SystemExit(str(exc))
    feature_types = set(args.feature_types) if args.feature_types else None
    _warn_if_index_unusable(
        args.input,
        args.index,
        feature_types=feature_types,
        attribute_filters=where_filters,
    )
    output = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        extract(
            source=args.input,
            fields=args.fields,
            feature_types=feature_types,
            index_path=args.index,
            attribute_filters=where_filters,
            output=output,
            separator=args.separator,
            deduplicate=not args.no_dedup,
        )
    finally:
        if output is not sys.stdout:
            output.close()


def _run_bed(args: argparse.Namespace) -> None:
    from .bed import to_bed

    try:
        where_filters = _parse_where_filters(args.where)
    except ValueError as exc:
        raise SystemExit(str(exc))
    _warn_if_index_unusable(
        args.input,
        args.index,
        feature_types={args.feature_type},
        attribute_filters=where_filters,
    )
    output = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        to_bed(
            source=args.input,
            feature_type=args.feature_type,
            id_fields=args.id_fields,
            id_separator=args.id_separator,
            index_path=args.index,
            attribute_filters=where_filters,
            output=output,
            deduplicate=not args.no_dedup,
        )
    finally:
        if output is not sys.stdout:
            output.close()


def _run_features(args: argparse.Namespace) -> None:
    from .parser import parse

    try:
        where_filters = _parse_where_filters(args.where)
    except ValueError as exc:
        raise SystemExit(str(exc))
    _warn_if_index_unusable(args.input, args.index, attribute_filters=where_filters)
    features: set[str] = set()
    for record in parse(args.input, index_path=args.index, attribute_filters=where_filters):
        features.add(record.feature)
    for f in sorted(features):
        print(f)


def _run_attributes(args: argparse.Namespace) -> None:
    from .parser import parse

    try:
        where_filters = _parse_where_filters(args.where)
    except ValueError as exc:
        raise SystemExit(str(exc))
    feature_types = set(args.feature_types) if args.feature_types else None
    _warn_if_index_unusable(
        args.input,
        args.index,
        feature_types=feature_types,
        attribute_filters=where_filters,
    )
    keys: set[str] = set()
    for record in parse(
        args.input,
        feature_types=feature_types,
        index_path=args.index,
        attribute_filters=where_filters,
    ):
        keys.update(record.attributes.keys())
    for k in sorted(keys):
        print(k)


def _run_index(args: argparse.Namespace) -> None:
    from .parser import build_index

    feature_types = set(args.feature_types) if args.feature_types else None
    attribute_keys = set(args.attribute_keys) if args.attribute_keys else None
    built = build_index(
        source=args.input,
        index_path=args.output,
        force=args.force,
        feature_types=feature_types,
        attribute_keys=attribute_keys,
    )
    print(built)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="gtf-parser",
        description="Parse GTF/GFF3 files: extract fields, convert to BED, and more.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subs = parser.add_subparsers(dest="command", required=True)

    # --- extract ---
    p_extract = subs.add_parser(
        "extract",
        help="Extract selected columns/attributes from a GTF/GFF3 file",
    )
    _add_common_args(p_extract)
    p_extract.add_argument(
        "-f",
        "--fields",
        nargs="+",
        required=True,
        help=(
            "Column names (seqname, source, feature, start, end, score, strand, frame) "
            "and/or attribute keys to extract"
        ),
    )
    p_extract.add_argument(
        "-s",
        "--separator",
        default="\t",
        help="Output column separator (default: tab)",
    )
    p_extract.set_defaults(func=_run_extract)

    # --- bed ---
    p_bed = subs.add_parser(
        "bed",
        help="Convert a GTF/GFF3 file to BED6 format",
    )
    p_bed.add_argument("input", help="Input GTF or GFF3 file (plain or .gz)")
    p_bed.add_argument(
        "feature_type",
        help="Feature type to convert (e.g. gene, exon, transcript)",
    )
    p_bed.add_argument(
        "-i",
        "--id-fields",
        nargs="+",
        required=True,
        help="Attribute keys / column names to combine into the BED name column",
    )
    p_bed.add_argument(
        "--id-separator",
        default="|",
        help="Separator for joining id-fields values (default: |)",
    )
    p_bed.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output file (default: stdout)",
    )
    p_bed.add_argument(
        "--no-dedup",
        action="store_true",
        help="Disable deduplication of output rows",
    )
    p_bed.add_argument(
        "--index",
        default=None,
        help="Optional SQLite index path created with `gtf-parser index`",
    )
    p_bed.add_argument(
        "--where",
        action="append",
        default=None,
        help=(
            "Attribute filter in the form key=value or key=v1,v2. "
            "Can be repeated to combine filters."
        ),
    )
    p_bed.set_defaults(func=_run_bed)

    # --- features ---
    p_features = subs.add_parser(
        "features",
        help="List all distinct feature types present in the file",
    )
    p_features.add_argument("input", help="Input GTF or GFF3 file")
    p_features.add_argument(
        "--index",
        default=None,
        help="Optional SQLite index path created with `gtf-parser index`",
    )
    p_features.add_argument(
        "--where",
        action="append",
        default=None,
        help=(
            "Attribute filter in the form key=value or key=v1,v2. "
            "Can be repeated to combine filters."
        ),
    )
    p_features.set_defaults(func=_run_features)

    # --- attributes ---
    p_attrs = subs.add_parser(
        "attributes",
        help="List all distinct attribute keys present in the file",
    )
    p_attrs.add_argument("input", help="Input GTF or GFF3 file")
    p_attrs.add_argument(
        "-t",
        "--feature-types",
        nargs="+",
        default=None,
        help="Restrict to these feature type(s)",
    )
    p_attrs.add_argument(
        "--index",
        default=None,
        help="Optional SQLite index path created with `gtf-parser index`",
    )
    p_attrs.add_argument(
        "--where",
        action="append",
        default=None,
        help=(
            "Attribute filter in the form key=value or key=v1,v2. "
            "Can be repeated to combine filters."
        ),
    )
    p_attrs.set_defaults(func=_run_attributes)

    # --- index ---
    p_index = subs.add_parser(
        "index",
        help="Build an SQLite index for faster repeated retrieval",
    )
    p_index.add_argument("input", help="Input GTF or GFF3 file (plain or .gz)")
    p_index.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output index path (default: <input>.idx.sqlite)",
    )
    p_index.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing index file if present",
    )
    p_index.add_argument(
        "-t",
        "--feature-types",
        nargs="+",
        default=None,
        help="Only index these feature type(s) for faster/smaller indexes",
    )
    p_index.add_argument(
        "-k",
        "--attribute-keys",
        nargs="+",
        default=None,
        help="Only index these attribute keys for --where lookups",
    )
    p_index.set_defaults(func=_run_index)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
