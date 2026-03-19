"""Tests for gtf_parser."""

import io
import textwrap

import pytest

from gtf_parser.parser import Record, build_index, index_status, parse
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

GTF_DATA_MULTI_GENE = textwrap.dedent("""\
    ##description: test
    chr1\thavana\tgene\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "ABC";
    chr1\thavana\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ABC";
    chr2\thavana\tgene\t10\t90\t.\t-\t.\tgene_id "G2"; gene_name "DEF";
    chr2\thavana\ttranscript\t10\t90\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "DEF";
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


def test_parse_attribute_filters_stream():
    records = list(
        parse(
            io.StringIO(GTF_DATA_MULTI_GENE),
            feature_types={"gene"},
            attribute_filters={"gene_id": {"G2"}},
        )
    )
    assert len(records) == 1
    assert records[0].get("gene_name") == "DEF"


def test_parse_attribute_filters_stream_multi_value():
    records = list(
        parse(
            io.StringIO(GTF_DATA_MULTI_GENE),
            feature_types={"gene"},
            attribute_filters={"gene_id": {"G1", "G2"}},
        )
    )
    assert len(records) == 2


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


def test_parse_with_index(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    records = list(parse(source, feature_types={"gene"}, index_path=index))
    assert len(records) == 1
    assert records[0].feature == "gene"
    assert records[0].get("gene_id") == "G1"


def test_index_status_ok(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    usable, reason = index_status(source, index)
    assert usable
    assert reason == "ok"


def test_index_status_feature_subset_missing(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source, feature_types={"gene"})

    usable, reason = index_status(source, index, feature_types={"exon"})
    assert not usable
    assert "missing feature type" in reason


def test_index_status_attribute_subset_missing(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source, attribute_keys={"gene_id"})

    usable, reason = index_status(
        source,
        index,
        attribute_filters={"transcript_id": {"T1"}},
    )
    assert not usable
    assert "missing attribute key" in reason


def test_parse_with_index_and_attribute_filters(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA_MULTI_GENE, encoding="utf-8")
    index = build_index(source)

    records = list(
        parse(
            source,
            feature_types={"gene"},
            index_path=index,
            attribute_filters={"gene_id": {"G2"}},
        )
    )
    assert len(records) == 1
    assert records[0].get("gene_name") == "DEF"


def test_parse_fallback_when_feature_not_indexed(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source, feature_types={"gene"})

    records = list(parse(source, feature_types={"exon"}, index_path=index))
    assert len(records) == 3
    assert all(record.feature == "exon" for record in records)


def test_extract_with_index(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    out = io.StringIO()
    extract(
        source,
        fields=["gene_id", "transcript_id"],
        feature_types={"exon"},
        output=out,
        index_path=index,
    )
    lines = out.getvalue().strip().split("\n")
    assert lines[0] == "gene_id\ttranscript_id"
    assert len(lines) == 3


def test_extract_with_index_and_attribute_filters(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA_MULTI_GENE, encoding="utf-8")
    index = build_index(source)

    out = io.StringIO()
    extract(
        source,
        fields=["gene_id", "gene_name"],
        feature_types={"gene"},
        output=out,
        index_path=index,
        attribute_filters={"gene_id": {"G2"}},
    )
    lines = out.getvalue().strip().split("\n")
    assert lines[0] == "gene_id\tgene_name"
    assert lines[1] == "G2\tDEF"
    assert len(lines) == 2


def test_bed_with_index(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    out = io.StringIO()
    to_bed(
        source,
        feature_type="gene",
        id_fields=["gene_id", "gene_name"],
        output=out,
        index_path=index,
    )
    lines = out.getvalue().strip().split("\n")
    assert len(lines) == 1
    assert lines[0].split("\t")[3] == "G1|ABC"


def test_bed_with_index_and_attribute_filters(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA_MULTI_GENE, encoding="utf-8")
    index = build_index(source)

    out = io.StringIO()
    to_bed(
        source,
        feature_type="gene",
        id_fields=["gene_id", "gene_name"],
        output=out,
        index_path=index,
        attribute_filters={"gene_id": {"G2"}},
    )
    lines = out.getvalue().strip().split("\n")
    assert len(lines) == 1
    assert lines[0].split("\t")[3] == "G2|DEF"


def test_index_stale_fallback(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    source.write_text(
        GTF_DATA + 'chr2\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G2"; gene_name "DEF";\n',
        encoding="utf-8",
    )

    records = list(parse(source, feature_types={"gene"}, index_path=index))
    assert len(records) == 2
    assert records[1].get("gene_id") == "G2"


def test_index_stale_fallback_with_attribute_filters(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    source.write_text(
        GTF_DATA + 'chr2\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G2"; gene_name "DEF";\n',
        encoding="utf-8",
    )

    records = list(
        parse(
            source,
            feature_types={"gene"},
            index_path=index,
            attribute_filters={"gene_id": {"G2"}},
        )
    )
    assert len(records) == 1
    assert records[0].get("gene_name") == "DEF"


def test_index_status_stale(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)
    source.write_text(GTF_DATA_MULTI_GENE, encoding="utf-8")

    usable, reason = index_status(source, index)
    assert not usable
    assert "stale" in reason


def test_build_index_force(tmp_path):
    source = tmp_path / "test.gtf"
    source.write_text(GTF_DATA, encoding="utf-8")
    index = build_index(source)

    with pytest.raises(FileExistsError):
        build_index(source, index_path=index)

    rebuilt = build_index(source, index_path=index, force=True)
    assert rebuilt == index
