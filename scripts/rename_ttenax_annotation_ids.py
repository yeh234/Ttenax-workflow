#!/usr/bin/env python3
"""Rename T. tenax annotation IDs using Ttenax_gene_id_mapping_v1.tsv.

Examples:
  python3 scripts/rename_ttenax_annotation_ids.py \
    --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
    --input old_annotation.gff3 \
    --output annotation/Ttenax_annotation_v1.gff3

  python3 scripts/rename_ttenax_annotation_ids.py \
    --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
    --input old_proteins.faa \
    --output annotation/Ttenax_predicted_proteins_v1.faa \
    --mode protein-fasta
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
        required = {
            "old_gene_id",
            "old_transcript_id",
            "old_protein_id",
            "new_gene_id",
            "new_transcript_id",
            "new_protein_id",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Mapping file is missing columns: {', '.join(sorted(missing))}")
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


def rename_text(
    input_path: Path,
    output_path: Path,
    gene_map: dict[str, str],
    transcript_map: dict[str, str],
    protein_map: dict[str, str],
) -> None:
    with input_path.open("r", encoding="utf-8", errors="replace", newline="") as inp, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as out:
        for line in inp:
            out.write(replace_ids(line, gene_map, transcript_map, protein_map))


def rename_fasta(input_path: Path, output_path: Path, id_map: dict[str, str]) -> None:
    with input_path.open("r", encoding="utf-8", errors="replace") as inp, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as out:
        for line in inp:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                first, sep, rest = header.partition(" ")
                out.write(f">{id_map.get(first, first)}{sep}{rest}\n")
            else:
                out.write(line)


def infer_mode(path: Path) -> str:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".faa") or suffixes.endswith(".aa.fasta"):
        return "protein-fasta"
    if suffixes.endswith(".fna") or suffixes.endswith(".nt.fasta"):
        return "cds-fasta"
    return "text"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=["auto", "text", "protein-fasta", "cds-fasta"],
        default="auto",
        help="Use text for GFF3/GTF/TSV-like files, protein-fasta for g1t1 headers, cds-fasta for g1.t1 headers.",
    )
    args = parser.parse_args()

    gene_map, transcript_map, protein_map = load_mapping(args.mapping)
    mode = infer_mode(args.input) if args.mode == "auto" else args.mode
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if mode == "protein-fasta":
        rename_fasta(args.input, args.output, protein_map)
    elif mode == "cds-fasta":
        rename_fasta(args.input, args.output, transcript_map)
    else:
        rename_text(args.input, args.output, gene_map, transcript_map, protein_map)


if __name__ == "__main__":
    main()
