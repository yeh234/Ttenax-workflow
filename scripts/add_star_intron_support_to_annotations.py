#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_gff3_attrs(text):
    attrs = {}
    if text == ".":
        return attrs
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def format_gff3_attrs(attrs):
    if not attrs:
        return "."
    return ";".join(f"{key}={value}" for key, value in attrs.items())


def parse_gtf_attrs(text):
    attrs = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
        attrs[key] = value
    return attrs


def append_gtf_attrs(text, additions):
    base = text.strip()
    if base == ".":
        base = ""
    if base and not base.endswith(";"):
        base += ";"
    for key, value in additions.items():
        base += f' {key} "{value}";'
    return base.strip()


def yesno(condition):
    return "yes" if condition else "no"


def gene_from_transcript(transcript_id):
    return re.sub(r"\.t\d+$", "", transcript_id)


def read_supported_introns(path):
    supported = {}
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if len(row) < 15:
                raise ValueError(f"Expected 15 columns in {path}, got {len(row)}: {row}")
            intron_id = row[6]
            supported[intron_id] = {
                "key": row[0],
                "seqid": row[1],
                "start": int(row[2]),
                "end": int(row[3]),
                "strand": row[4],
                "transcript_id": row[5],
                "intron_id": intron_id,
                "star_unique_reads": int(row[11]),
                "star_multimap_reads": int(row[12]),
                "star_total_reads": int(row[13]),
                "star_annotated_junction": row[14],
            }
    return supported


def collect_gff3_models(gff3_path):
    transcript_to_gene = {}
    introns_by_transcript = defaultdict(list)
    intron_records = {}

    with open(gff3_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            attrs = parse_gff3_attrs(parts[8])
            if feature == "mRNA":
                transcript_id = attrs.get("ID")
                gene_id = attrs.get("Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
            elif feature == "intron":
                intron_id = attrs.get("ID")
                transcript_id = attrs.get("Parent")
                if not intron_id or not transcript_id:
                    continue
                record = {
                    "seqid": parts[0],
                    "source": parts[1],
                    "feature": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attrs": attrs,
                    "intron_id": intron_id,
                    "transcript_id": transcript_id,
                }
                introns_by_transcript[transcript_id].append(record)
                intron_records[intron_id] = record

    for transcript_id in introns_by_transcript:
        transcript_to_gene.setdefault(transcript_id, gene_from_transcript(transcript_id))

    return transcript_to_gene, introns_by_transcript, intron_records


def build_summaries(transcript_to_gene, introns_by_transcript, supported, supported_min5):
    transcript_summary = {}
    gene_summary = defaultdict(lambda: {
        "annotated_introns": 0,
        "star_supported_introns": 0,
        "star_supported_introns_min5": 0,
        "transcripts_with_introns": 0,
        "transcripts_with_star_intron_support": 0,
        "transcripts_all_introns_supported": 0,
    })

    for transcript_id, introns in introns_by_transcript.items():
        total = len(introns)
        support = sum(1 for intron in introns if intron["intron_id"] in supported)
        support_min5 = sum(1 for intron in introns if intron["intron_id"] in supported_min5)
        all_supported = total > 0 and support == total
        transcript_summary[transcript_id] = {
            "annotated_introns": total,
            "star_supported_introns": support,
            "star_supported_introns_min5": support_min5,
            "star_any_intron_support": yesno(support > 0),
            "star_all_introns_supported": yesno(all_supported),
        }

        gene_id = transcript_to_gene.get(transcript_id, gene_from_transcript(transcript_id))
        gs = gene_summary[gene_id]
        gs["annotated_introns"] += total
        gs["star_supported_introns"] += support
        gs["star_supported_introns_min5"] += support_min5
        gs["transcripts_with_introns"] += 1
        gs["transcripts_with_star_intron_support"] += int(support > 0)
        gs["transcripts_all_introns_supported"] += int(all_supported)

    for gene_id, summary in gene_summary.items():
        total = summary["annotated_introns"]
        support = summary["star_supported_introns"]
        summary["star_any_intron_support"] = yesno(support > 0)
        summary["star_all_introns_supported"] = yesno(total > 0 and support == total)

    return transcript_summary, dict(gene_summary)


def star_intron_attrs(intron_id, supported, supported_min5):
    support = supported.get(intron_id)
    if support is None:
        return {
            "STAR_junction_support": "no",
            "STAR_unique_reads": "0",
            "STAR_multimap_reads": "0",
            "STAR_total_reads": "0",
            "STAR_annotated_junction": "0",
            "STAR_min5_support": "no",
        }
    return {
        "STAR_junction_support": "exact",
        "STAR_unique_reads": str(support["star_unique_reads"]),
        "STAR_multimap_reads": str(support["star_multimap_reads"]),
        "STAR_total_reads": str(support["star_total_reads"]),
        "STAR_annotated_junction": str(support["star_annotated_junction"]),
        "STAR_min5_support": yesno(intron_id in supported_min5),
    }


def star_model_attrs(summary):
    return {
        "STAR_annotated_introns": str(summary["annotated_introns"]),
        "STAR_supported_introns": str(summary["star_supported_introns"]),
        "STAR_supported_introns_min5": str(summary["star_supported_introns_min5"]),
        "STAR_any_intron_support": summary["star_any_intron_support"],
        "STAR_all_introns_supported": summary["star_all_introns_supported"],
    }


def write_gff3_with_support(in_path, out_path, transcript_summary, gene_summary,
                            supported, supported_min5):
    with open(in_path, "r", encoding="utf-8", errors="replace") as inp, \
            open(out_path, "w", encoding="utf-8", newline="") as out:
        for line in inp:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                out.write(line)
                continue
            feature = parts[2]
            attrs = parse_gff3_attrs(parts[8])
            if feature == "gene":
                gene_id = attrs.get("ID")
                if gene_id in gene_summary:
                    attrs.update(star_model_attrs(gene_summary[gene_id]))
            elif feature == "mRNA":
                transcript_id = attrs.get("ID")
                if transcript_id in transcript_summary:
                    attrs.update(star_model_attrs(transcript_summary[transcript_id]))
            elif feature == "intron":
                intron_id = attrs.get("ID")
                if intron_id:
                    attrs.update(star_intron_attrs(intron_id, supported, supported_min5))
            parts[8] = format_gff3_attrs(attrs)
            out.write("\t".join(parts) + "\n")


def gtf_intron_attrs(record, transcript_to_gene, supported, supported_min5):
    transcript_id = record["transcript_id"]
    gene_id = transcript_to_gene.get(transcript_id, gene_from_transcript(transcript_id))
    additions = {
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "intron_id": record["intron_id"],
    }
    for key, value in star_intron_attrs(record["intron_id"], supported, supported_min5).items():
        additions[key.lower()] = value
    return " ".join(f'{key} "{value}";' for key, value in additions.items())


def write_gtf_with_support(in_path, out_path, transcript_to_gene, introns_by_transcript,
                           transcript_summary, supported, supported_min5):
    with open(in_path, "r", encoding="utf-8", errors="replace") as inp, \
            open(out_path, "w", encoding="utf-8", newline="") as out:
        for line in inp:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                out.write(line)
                continue
            attrs = parse_gtf_attrs(parts[8])
            transcript_id = attrs.get("transcript_id")
            if parts[2] == "transcript" and transcript_id in transcript_summary:
                additions = {
                    key.lower(): value
                    for key, value in star_model_attrs(transcript_summary[transcript_id]).items()
                }
                parts[8] = append_gtf_attrs(parts[8], additions)
                out.write("\t".join(parts) + "\n")

                supported_introns = [
                    intron
                    for intron in introns_by_transcript[transcript_id]
                    if intron["intron_id"] in supported
                ]
                for intron in sorted(supported_introns,
                                     key=lambda x: (x["start"], x["end"], x["intron_id"])):
                    intron_parts = [
                        intron["seqid"],
                        intron["source"],
                        "intron",
                        str(intron["start"]),
                        str(intron["end"]),
                        intron["score"],
                        intron["strand"],
                        ".",
                        gtf_intron_attrs(intron, transcript_to_gene, supported, supported_min5),
                    ]
                    out.write("\t".join(intron_parts) + "\n")
            else:
                out.write(line)


def write_supported_introns_table(out_path, transcript_to_gene, supported, supported_min5):
    fieldnames = [
        "gene_id",
        "transcript_id",
        "intron_id",
        "seqid",
        "start",
        "end",
        "strand",
        "star_unique_reads",
        "star_multimap_reads",
        "star_total_reads",
        "star_annotated_junction",
        "star_min5_support",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for intron_id in sorted(supported):
            item = supported[intron_id]
            transcript_id = item["transcript_id"]
            writer.writerow({
                "gene_id": transcript_to_gene.get(transcript_id, gene_from_transcript(transcript_id)),
                "transcript_id": transcript_id,
                "intron_id": intron_id,
                "seqid": item["seqid"],
                "start": item["start"],
                "end": item["end"],
                "strand": item["strand"],
                "star_unique_reads": item["star_unique_reads"],
                "star_multimap_reads": item["star_multimap_reads"],
                "star_total_reads": item["star_total_reads"],
                "star_annotated_junction": item["star_annotated_junction"],
                "star_min5_support": yesno(intron_id in supported_min5),
            })


def write_summary(out_path, transcript_summary, gene_summary, intron_records,
                  supported, supported_min5, gff3_out, gtf_out, supported_table):
    rows = [
        ("annotated_introns_in_gff3", len(intron_records), "Input GFF3 intron features"),
        ("star_exact_supported_introns", len(supported), "Exact coordinate+strand STAR junction matches"),
        ("star_exact_supported_introns_min5", len(supported_min5), "Exact STAR matches with total junction reads >=5"),
        ("explicit_introns_written_to_gtf", len(supported), "Only STAR-supported intron features are written to GTF"),
        ("intron_containing_transcripts", len(transcript_summary), "Transcripts with >=1 annotated intron"),
        ("transcripts_with_star_supported_introns",
         sum(1 for s in transcript_summary.values() if s["star_supported_introns"] > 0),
         "Transcripts with >=1 STAR-supported intron"),
        ("transcripts_all_introns_supported",
         sum(1 for s in transcript_summary.values() if s["star_all_introns_supported"] == "yes"),
         "Transcripts whose annotated introns are all STAR-supported"),
        ("intron_containing_genes", len(gene_summary), "Genes with >=1 annotated intron"),
        ("genes_with_star_supported_introns",
         sum(1 for s in gene_summary.values() if s["star_supported_introns"] > 0),
         "Genes with >=1 STAR-supported intron"),
        ("genes_all_introns_supported",
         sum(1 for s in gene_summary.values() if s["star_all_introns_supported"] == "yes"),
         "Genes whose annotated introns are all STAR-supported"),
        ("output_gff3", gff3_out, "GFF3 with STAR support attributes"),
        ("output_gtf", gtf_out, "GTF with transcript STAR support attributes and STAR-supported intron rows only"),
        ("output_supported_introns_table", supported_table, "Supported intron table with gene IDs"),
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["metric", "value", "note"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Add exact STAR splice-junction support annotations to T. tenax GFF3/GTF files."
    )
    parser.add_argument("--gff3", default="Express/Tt_Final_Annotation.gff3")
    parser.add_argument("--gtf", default="Express/Tt_Final_Annotation.gtf")
    parser.add_argument("--supported", default="Express/intron_check/annotation_introns_supported_by_STAR.tsv")
    parser.add_argument("--supported-min5", default="Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv")
    parser.add_argument("--out-prefix", default="Express/Tt_Final_Annotation.STAR_intron_support")
    args = parser.parse_args()

    gff3_path = Path(args.gff3)
    gtf_path = Path(args.gtf)
    supported_path = Path(args.supported)
    supported_min5_path = Path(args.supported_min5)
    for path in [gff3_path, gtf_path, supported_path, supported_min5_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    supported = read_supported_introns(supported_path)
    supported_min5 = read_supported_introns(supported_min5_path)
    transcript_to_gene, introns_by_transcript, intron_records = collect_gff3_models(gff3_path)
    transcript_summary, gene_summary = build_summaries(
        transcript_to_gene, introns_by_transcript, supported, supported_min5
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_gff3 = Path(str(out_prefix) + ".gff3")
    out_gtf = Path(str(out_prefix) + ".gtf")
    out_summary = Path(str(out_prefix) + ".summary.tsv")
    out_supported_table = Path(str(out_prefix) + ".supported_introns.tsv")

    write_gff3_with_support(gff3_path, out_gff3, transcript_summary, gene_summary,
                            supported, supported_min5)
    write_gtf_with_support(gtf_path, out_gtf, transcript_to_gene, introns_by_transcript,
                           transcript_summary, supported, supported_min5)
    write_supported_introns_table(out_supported_table, transcript_to_gene,
                                  supported, supported_min5)
    write_summary(out_summary, transcript_summary, gene_summary, intron_records,
                  supported, supported_min5, out_gff3, out_gtf, out_supported_table)

    print(f"Wrote: {out_gff3}")
    print(f"Wrote: {out_gtf}")
    print(f"Wrote: {out_summary}")
    print(f"Wrote: {out_supported_table}")


if __name__ == "__main__":
    main()
