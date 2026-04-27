#!/usr/bin/env python3
"""Apply T. tenax v1 gene ID mapping to tabular files.

By default, every cell is scanned for old gene, transcript, protein, and
subfeature IDs. Use --columns to restrict replacement to selected columns.

Example:
  python3 scripts/apply_gene_id_mapping_to_tables.py \
    --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
    --input old_table.tsv \
    --output new_table.tsv \
    --add-old-gene-id
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ID_RE = re.compile(
    r"\bg\d+\.t\d+\.(?:CDS|exon|intron|start|stop)\d+\b"
    r"|\bg\d+\.t\d+\b"
    r"|\bg\d+t\d+\b"
    r"|\bg\d+\b"
)


def load_mapping(path: Path) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    gene_map: dict[str, str] = {}
    transcript_map: dict[str, str] = {}
    protein_map: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_map[row["old_gene_id"]] = row["new_gene_id"]
            transcript_map[row["old_transcript_id"]] = row["new_transcript_id"]
            protein_map[row["old_protein_id"]] = row["new_protein_id"]
    return gene_map, transcript_map, protein_map


def replace_ids(
    text: str,
    gene_map: dict[str, str],
    transcript_map: dict[str, str],
    protein_map: dict[str, str],
) -> str:
    def convert(match: re.Match[str]) -> str:
        token = match.group(0)
        subfeature = re.match(r"^(g\d+\.t\d+)\.((?:CDS|exon|intron|start|stop)\d+)$", token)
        if subfeature:
            transcript_id, suffix = subfeature.groups()
            return f"{transcript_map.get(transcript_id, transcript_id)}.{suffix}"
        if token in transcript_map:
            return transcript_map[token]
        if token in protein_map:
            return protein_map[token]
        if token in gene_map:
            return gene_map[token]
        return token

    return ID_RE.sub(convert, text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--columns",
        help="Comma-separated column names to transform. Default: transform all columns.",
    )
    parser.add_argument(
        "--add-old-gene-id",
        action="store_true",
        help="If GeneID is mapped to AC16WH_, add an old_gene_id traceability column after GeneID.",
    )
    args = parser.parse_args()

    gene_map, transcript_map, protein_map = load_mapping(args.mapping)
    reverse_gene_map = {new: old for old, new in gene_map.items()}
    selected = set(args.columns.split(",")) if args.columns else None

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(newline="") as inp, args.output.open("w", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Input table has no header")

        fieldnames = list(reader.fieldnames)
        if args.add_old_gene_id and "GeneID" in fieldnames and "old_gene_id" not in fieldnames:
            gene_index = fieldnames.index("GeneID")
            fieldnames.insert(gene_index + 1, "old_gene_id")

        writer = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            original_gene = row.get("GeneID", "")
            for field in reader.fieldnames:
                if selected is None or field in selected:
                    row[field] = replace_ids(row.get(field, ""), gene_map, transcript_map, protein_map)
            if args.add_old_gene_id and "GeneID" in row:
                row["old_gene_id"] = reverse_gene_map.get(row["GeneID"], original_gene)
            writer.writerow(row)


if __name__ == "__main__":
    main()
