# gtf-parser

Parse GTF and GFF3 files — extract arbitrary fields, deduplicate, and convert
to BED format. Pure Python, no external dependencies.

## Installation

```bash
pip install -e .
```

## Quick start

### List available feature types

```bash
gtf-parser features genes.gtf
```

### List attribute keys

```bash
gtf-parser attributes genes.gtf
gtf-parser attributes genes.gtf -t gene   # restrict to "gene" features
```

### Extract fields

Extract any combination of main columns (`seqname`, `source`, `feature`,
`start`, `end`, `score`, `strand`, `frame`) and annotation attribute keys.
Output is deduplicated by default.

```bash
# Extract gene_id and gene_name for all "gene" features
gtf-parser extract genes.gtf -t gene -f gene_id gene_name

# Extract from all feature types, include coordinates
gtf-parser extract genes.gtf -f seqname start end gene_id transcript_id

# Write to a file instead of stdout
gtf-parser extract genes.gtf -t exon -f gene_id transcript_id exon_number -o exons.tsv

# Keep duplicates
gtf-parser extract genes.gtf -t exon -f gene_id --no-dedup
```

### Convert to BED

Convert features of a given type to BED6 format. The `name` column is built
from a configurable combination of attribute keys / columns.

```bash
# BED with name = gene_id|gene_name
gtf-parser bed genes.gtf gene -i gene_id gene_name

# Custom separator
gtf-parser bed genes.gtf exon -i gene_id transcript_id exon_number --id-separator ":"

# Write to file
gtf-parser bed genes.gtf transcript -i transcript_id -o transcripts.bed
```

## Supported formats

- **GTF** (Gene Transfer Format) — attributes like `gene_id "ENSG00000...";`
- **GFF3** (General Feature Format v3) — attributes like `ID=gene0001;Name=...`
- Gzip-compressed files (`.gtf.gz`, `.gff3.gz`) are handled transparently.

The format is auto-detected from the first data line.

## Python API

```python
from gtf_parser.parser import parse

for record in parse("genes.gtf", feature_types={"gene"}):
    print(record.seqname, record.get("gene_id"), record.get("gene_name"))
```

```python
from gtf_parser.extract import extract

extract("genes.gtf", fields=["gene_id", "gene_name"], feature_types={"gene"})
```

```python
from gtf_parser.bed import to_bed

to_bed("genes.gtf", feature_type="gene", id_fields=["gene_id", "gene_name"])
```

## Contributors

- [Edoardo Giacopuzzi](https://github.com/edg1983)

With help from Claude Opus 4.6.