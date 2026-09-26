# Release Naming Notes

This repository uses the *T. tenax* annotation v1.0.0 release naming scheme.

## Prefix

All public release gene IDs use:

```text
AC16WH_
```

Examples:

```text
Gene:       AC16WH_000001
Transcript: AC16WH_000001.t1
Protein:    AC16WH_000001.p1
```

## Mapping

The file `annotation/Ttenax_gene_id_mapping_v1.tsv` maps the original GALBA-style IDs to release IDs:

```text
old_gene_id  old_transcript_id  old_protein_id  new_gene_id     new_transcript_id  new_protein_id
g1           g1.t1              g1t1            AC16WH_000001   AC16WH_000001.t1  AC16WH_000001.p1
```

The frozen release matrix `annotation/Ttenax_multiomics_gene_evidence_matrix_v1.tsv` keeps:

- `GeneID`: release ID
- `old_gene_id`: original ID for traceability

## Files Using Release IDs

```text
annotation/Ttenax_annotation_v1.gff3
annotation/Ttenax_annotation_v1.gtf
annotation/Ttenax_predicted_proteins_v1.faa
annotation/Ttenax_predicted_cds_v1.fna
annotation/Ttenax_multiomics_gene_evidence_matrix_v1.tsv
annotation/star_intron_support/Ttenax_STAR_supported_introns_v1.tsv
```

The full STAR-support-enhanced GFF3/GTF files use the same release IDs but are
kept in `Ttenax_annotation_release_v1.0.0.tar.gz` rather than duplicated in this
compact GitHub repository.

## Scripts

- `scripts/rename_ttenax_annotation_ids.py` renames GFF3/GTF/FASTA annotation IDs.
- `scripts/apply_gene_id_mapping_to_tables.py` applies the same mapping to TSV-style tables.
- `scripts/build_multiomics_gene_evidence_matrix.py` rebuilds the evidence matrix from featureCounts and proteomics inputs.

## Versioning

Use the `v1` suffix for files inside the annotation release. The archived GitHub/Zenodo release remains `v1.0.0`.

The September 2026 revision on `main` uses plotting-script version `20260916.1` and the updated workbook name `Supplementary_Tables_S1-S7_20260926.xlsx`. Its exported figure metrics retain the dated `Figure*_20260916_v1_*.tsv` filenames because their contents and plotting version are unchanged. Do not interpret the script version or these filename suffixes as a newly published Zenodo release. No v1.1 release tag has been created for this synchronization.
