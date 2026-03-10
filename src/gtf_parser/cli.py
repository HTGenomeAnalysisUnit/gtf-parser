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


def _run_extract(args: argparse.Namespace) -> None:
    from .extract import extract

    output = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        feature_types = set(args.feature_types) if args.feature_types else None
        extract(
            source=args.input,
            fields=args.fields,
            feature_types=feature_types,
            output=output,
            separator=args.separator,
            deduplicate=not args.no_dedup,
        )
    finally:
        if output is not sys.stdout:
            output.close()


def _run_bed(args: argparse.Namespace) -> None:
    from .bed import to_bed

    output = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        to_bed(
            source=args.input,
            feature_type=args.feature_type,
            id_fields=args.id_fields,
            id_separator=args.id_separator,
            output=output,
            deduplicate=not args.no_dedup,
        )
    finally:
        if output is not sys.stdout:
            output.close()


def _run_features(args: argparse.Namespace) -> None:
    from .parser import parse

    features: set[str] = set()
    for record in parse(args.input):
        features.add(record.feature)
    for f in sorted(features):
        print(f)


def _run_attributes(args: argparse.Namespace) -> None:
    from .parser import parse

    feature_types = set(args.feature_types) if args.feature_types else None
    keys: set[str] = set()
    for record in parse(args.input, feature_types=feature_types):
        keys.update(record.attributes.keys())
    for k in sorted(keys):
        print(k)


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
    p_bed.set_defaults(func=_run_bed)

    # --- features ---
    p_features = subs.add_parser(
        "features",
        help="List all distinct feature types present in the file",
    )
    p_features.add_argument("input", help="Input GTF or GFF3 file")
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
    p_attrs.set_defaults(func=_run_attributes)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
