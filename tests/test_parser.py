"""Tests for gtf_parser."""

import io
import textwrap

from gtf_parser.parser import Record, parse
from gtf_parser.extract import extract
from gtf_parser.bed import to_bed

GTF_DATA = textwrap.dedent("""\
    ##description: test
    chr1\thavana\tgene\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "ABC";
    chr1\thavana\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ABC";
    chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
    chr1\thavana\texon\t300\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";
    chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
""")

GFF3_DATA = textwrap.dedent("""\
    ##gff-version 3
    chr1\tensembl\tgene\t100\t500\t.\t+\t.\tID=gene0001;Name=ABC
    chr1\tensembl\tmRNA\t100\t500\t.\t+\t.\tID=mRNA0001;Parent=gene0001
""")


def test_parse_gtf_all():
    records = list(parse(io.StringIO(GTF_DATA)))
    assert len(records) == 5
    assert records[0].feature == "gene"
    assert records[0].get("gene_id") == "G1"


def test_parse_gtf_filter():
    records = list(parse(io.StringIO(GTF_DATA), feature_types={"exon"}))
    assert len(records) == 3
    assert all(r.feature == "exon" for r in records)


def test_parse_gff3():
    records = list(parse(io.StringIO(GFF3_DATA)))
    assert len(records) == 2
    assert records[0].get("ID") == "gene0001"
    assert records[0].get("Name") == "ABC"


def test_extract_dedup():
    out = io.StringIO()
    extract(io.StringIO(GTF_DATA), fields=["gene_id", "transcript_id"], feature_types={"exon"}, output=out)
    lines = out.getvalue().strip().split("\n")
    # header + 2 unique rows (the duplicate exon_number=1 row is deduped)
    assert lines[0] == "gene_id\ttranscript_id"
    assert len(lines) == 3  # header + 2 data rows


def test_extract_no_dedup():
    out = io.StringIO()
    extract(io.StringIO(GTF_DATA), fields=["gene_id"], feature_types={"exon"}, output=out, deduplicate=False)
    lines = out.getvalue().strip().split("\n")
    assert len(lines) == 4  # header + 3 data rows


def test_bed_conversion():
    out = io.StringIO()
    to_bed(io.StringIO(GTF_DATA), feature_type="gene", id_fields=["gene_id", "gene_name"], output=out)
    lines = out.getvalue().strip().split("\n")
    assert len(lines) == 1
    parts = lines[0].split("\t")
    assert parts[0] == "chr1"
    assert parts[1] == "99"  # 0-based
    assert parts[2] == "500"
    assert parts[3] == "G1|ABC"
    assert parts[5] == "+"


def test_bed_dedup():
    out = io.StringIO()
    to_bed(io.StringIO(GTF_DATA), feature_type="exon", id_fields=["gene_id", "exon_number"], output=out)
    lines = out.getvalue().strip().split("\n")
    # 3 exon records but one is a duplicate → 2 unique BED lines
    assert len(lines) == 2


def test_record_get_column():
    r = Record(seqname="chr1", source="test", feature="gene", start=1, end=100, score=".", strand="+", frame=".")
    assert r.get("seqname") == "chr1"
    assert r.get("start") == "1"
    assert r.get("nonexistent") is None
