#!/usr/bin/env python3
"""Build the T. tenax multi-omics gene evidence matrix.

Inputs are unstranded featureCounts gene tables for DRS, short-read cDNA, and
long-read cDNA, plus an optional gene-level proteomics table. The output uses
AC16WH_ release gene IDs and retains the original GALBA-style gene ID in
old_gene_id. Proteomics best-protein IDs are left as reported by default to
match the frozen v1 release table; use --map-protein-ids to convert those to
AC16WH_ protein IDs.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path


OUTPUT_COLUMNS = [
    "GeneID",
    "old_gene_id",
    "DRS_count",
    "DRS_support",
    "DRS_ge5",
    "DRS_ge10",
    "SR_cDNA_count",
    "SR_cDNA_support",
    "SR_cDNA_ge5",
    "SR_cDNA_ge10",
    "LR_cDNA_count",
    "LR_cDNA_support",
    "LR_cDNA_ge5",
    "LR_cDNA_ge10",
    "RNA_platforms_supported",
    "Any_RNA_support",
    "All_three_RNA_support",
    "Proteomics_support",
    "Max_unique_peptides",
    "Sum_unique_peptides_reported",
    "Proteomics_total_PSMs",
    "Proteomics_total_peptides_reported",
    "Proteomics_max_coverage_percent",
    "Proteomics_max_Mascot_score",
    "Best_protein_by_unique_peptides",
    "High_confidence_proteomics_2_unique_peptides",
    "Strong_proteomics_3_unique_peptides",
    "Very_strong_proteomics_5_unique_peptides",
    "Multiomics_support",
    "Evidence_class",
]


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def load_mapping(path: Path) -> tuple[OrderedDict[str, str], dict[str, str]]:
    gene_map: OrderedDict[str, str] = OrderedDict()
    protein_map: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_map.setdefault(row["old_gene_id"], row["new_gene_id"])
            protein_map[row["old_protein_id"]] = row["new_protein_id"]
    return gene_map, protein_map


def read_featurecounts(path: Path) -> OrderedDict[str, int]:
    counts: OrderedDict[str, int] = OrderedDict()
    with path.open(newline="") as handle:
        data_lines = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        if reader.fieldnames is None or "Geneid" not in reader.fieldnames:
            raise SystemExit(f"{path} does not look like a featureCounts table")
        count_column = reader.fieldnames[-1]
        for row in reader:
            raw = row[count_column]
            counts[row["Geneid"]] = int(float(raw)) if raw else 0
    return counts


def convert_protein_ids(text: str, protein_map: dict[str, str]) -> str:
    def convert(match: re.Match[str]) -> str:
        token = match.group(0)
        return protein_map.get(token, token)

    return re.sub(r"\bg\d+t\d+\b", convert, text or "")


def read_proteomics(
    path: Path | None,
    protein_map: dict[str, str],
    map_protein_ids: bool = False,
) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    proteomics: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row.get("GeneID", "")
            if not gene:
                continue
            proteomics[gene] = {
                "Proteomics_support": row.get("Proteomics_support", "yes") or "yes",
                "Max_unique_peptides": row.get("Max_unique_peptides", "0") or "0",
                "Sum_unique_peptides_reported": row.get("Sum_unique_peptides_reported", "0") or "0",
                "Proteomics_total_PSMs": row.get("Total_PSMs", row.get("Proteomics_total_PSMs", "0")) or "0",
                "Proteomics_total_peptides_reported": row.get(
                    "Total_peptides_reported", row.get("Proteomics_total_peptides_reported", "0")
                )
                or "0",
                "Proteomics_max_coverage_percent": row.get("Max_coverage_percent", "0") or "0",
                "Proteomics_max_Mascot_score": row.get("Max_Mascot_score", "0") or "0",
                "Best_protein_by_unique_peptides": (
                    convert_protein_ids(row.get("Best_protein_by_unique_peptides", ""), protein_map)
                    if map_protein_ids
                    else row.get("Best_protein_by_unique_peptides", "")
                ),
                "High_confidence_proteomics_2_unique_peptides": row.get(
                    "High_confidence_2_unique_peptides",
                    row.get("High_confidence_proteomics_2_unique_peptides", "no"),
                )
                or "no",
                "Strong_proteomics_3_unique_peptides": row.get(
                    "Strong_support_3_unique_peptides",
                    row.get("Strong_proteomics_3_unique_peptides", "no"),
                )
                or "no",
                "Very_strong_proteomics_5_unique_peptides": row.get(
                    "Very_strong_support_5_unique_peptides",
                    row.get("Very_strong_proteomics_5_unique_peptides", "no"),
                )
                or "no",
            }
    return proteomics


def evidence_class(rna_platforms: int, proteomics: bool, all_three: bool, high_conf: bool) -> str:
    if proteomics and all_three and high_conf:
        return "three_RNA_and_high_confidence_proteomics_supported"
    if proteomics and all_three:
        return "three_RNA_and_proteomics_supported"
    if proteomics:
        return "RNA_and_proteomics_supported"
    if all_three:
        return "three_RNA_supported"
    if rna_platforms >= 2:
        return "multi_RNA_supported"
    if rna_platforms == 1:
        return "single_RNA_supported"
    return "no_detected_omics_support"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--drs-counts", type=Path, required=True)
    parser.add_argument("--sr-cdna-counts", type=Path, required=True)
    parser.add_argument("--lr-cdna-counts", type=Path, required=True)
    parser.add_argument("--proteomics-gene-level", type=Path)
    parser.add_argument(
        "--map-protein-ids",
        action="store_true",
        help="Convert proteomics best-protein IDs from g1t1 style to AC16WH_000001.p1 style.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    gene_map, protein_map = load_mapping(args.mapping)
    drs = read_featurecounts(args.drs_counts)
    sr = read_featurecounts(args.sr_cdna_counts)
    lr = read_featurecounts(args.lr_cdna_counts)
    proteomics = read_proteomics(args.proteomics_gene_level, protein_map, args.map_protein_ids)

    old_gene_order = list(drs.keys())
    for gene in gene_map:
        if gene not in drs:
            old_gene_order.append(gene)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for old_gene in old_gene_order:
            if old_gene not in gene_map:
                continue
            drs_count = drs.get(old_gene, 0)
            sr_count = sr.get(old_gene, 0)
            lr_count = lr.get(old_gene, 0)
            platform_flags = [drs_count > 0, sr_count > 0, lr_count > 0]
            rna_platforms = sum(platform_flags)
            all_three = rna_platforms == 3

            prot = proteomics.get(old_gene, {})
            proteomics_support = prot.get("Proteomics_support", "no") == "yes"
            high_conf = prot.get("High_confidence_proteomics_2_unique_peptides", "no") == "yes"

            row = {
                "GeneID": gene_map[old_gene],
                "old_gene_id": old_gene,
                "DRS_count": drs_count,
                "DRS_support": yesno(drs_count > 0),
                "DRS_ge5": yesno(drs_count >= 5),
                "DRS_ge10": yesno(drs_count >= 10),
                "SR_cDNA_count": sr_count,
                "SR_cDNA_support": yesno(sr_count > 0),
                "SR_cDNA_ge5": yesno(sr_count >= 5),
                "SR_cDNA_ge10": yesno(sr_count >= 10),
                "LR_cDNA_count": lr_count,
                "LR_cDNA_support": yesno(lr_count > 0),
                "LR_cDNA_ge5": yesno(lr_count >= 5),
                "LR_cDNA_ge10": yesno(lr_count >= 10),
                "RNA_platforms_supported": rna_platforms,
                "Any_RNA_support": yesno(rna_platforms > 0),
                "All_three_RNA_support": yesno(all_three),
                "Proteomics_support": yesno(proteomics_support),
                "Max_unique_peptides": prot.get("Max_unique_peptides", "0"),
                "Sum_unique_peptides_reported": prot.get("Sum_unique_peptides_reported", "0"),
                "Proteomics_total_PSMs": prot.get("Proteomics_total_PSMs", "0"),
                "Proteomics_total_peptides_reported": prot.get("Proteomics_total_peptides_reported", "0"),
                "Proteomics_max_coverage_percent": prot.get("Proteomics_max_coverage_percent", "0"),
                "Proteomics_max_Mascot_score": prot.get("Proteomics_max_Mascot_score", "0"),
                "Best_protein_by_unique_peptides": prot.get("Best_protein_by_unique_peptides", ""),
                "High_confidence_proteomics_2_unique_peptides": prot.get(
                    "High_confidence_proteomics_2_unique_peptides", "no"
                ),
                "Strong_proteomics_3_unique_peptides": prot.get("Strong_proteomics_3_unique_peptides", "no"),
                "Very_strong_proteomics_5_unique_peptides": prot.get(
                    "Very_strong_proteomics_5_unique_peptides", "no"
                ),
                "Multiomics_support": yesno((rna_platforms > 0) and proteomics_support),
                "Evidence_class": evidence_class(rna_platforms, proteomics_support, all_three, high_conf),
            }
            writer.writerow(row)


if __name__ == "__main__":
    main()
