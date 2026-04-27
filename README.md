# Ttenax-workflow

Compact workflow and release-support repository for the *Trichomonas tenax* genome annotation v1.0.0.

This repository intentionally keeps only the v1 annotation release files, command-level SOPs, and small scripts needed to reproduce the annotation naming, multi-omics evidence matrix, and manuscript figure source workflows. Large raw or intermediate data such as POD5/FASTQ/BAM, STAR indexes, Flye directories, BUSCO directories, GALBA output, EggNOG output, and InterProScan output are not stored here.

## Layout

```text
annotation/
  Ttenax_annotation_v1.gff3
  Ttenax_annotation_v1.gtf
  Ttenax_predicted_proteins_v1.faa
  Ttenax_predicted_cds_v1.fna
  Ttenax_gene_id_mapping_v1.tsv
  Ttenax_multiomics_gene_evidence_matrix_v1.tsv
  star_intron_support/

scripts/
  rename_ttenax_annotation_ids.py
  apply_gene_id_mapping_to_tables.py
  build_multiomics_gene_evidence_matrix.py
  figure_plotting_scripts/

SOP/
  assembly_annotation_SOP.md
  multiomics_validation_SOP.md
  STAR_intron_support_SOP.md

tables/
  Supplementary_Tables_S1-S5.xlsx
  Figure1_source_data.tsv
  ...
  Figure6_source_data.tsv
```

## Release Files

The annotation files were imported from `Ttenax_annotation_release_v1.0.0.tar.gz`.

- Gene IDs use the `AC16WH_` prefix.
- Transcript IDs use `.t1`, for example `AC16WH_000001.t1`.
- Protein IDs use `.p1`, for example `AC16WH_000001.p1`.
- The original GALBA-style IDs are retained in `Ttenax_gene_id_mapping_v1.tsv` and in the `old_gene_id` column of the multi-omics evidence matrix.

## Main Workflows

The SOPs document the command-level workflow from Nanopore genomic FASTQ/BAM through Flye assembly, haplotig/depth diagnosis, Dorado polishing, BUSCO/MetaEuk assessment, GALBA structural annotation, functional annotation merging, and final reference/annotation release.

The multi-omics SOP documents DRS, long-read cDNA, short-read cDNA, and proteomics evidence generation. The included `build_multiomics_gene_evidence_matrix.py` script rebuilds the gene-level matrix from featureCounts outputs and the proteomics gene-level table.

The compact STAR support evidence table is included at `annotation/star_intron_support/Ttenax_STAR_supported_introns_v1.tsv`. The full STAR-support-enhanced GFF3/GTF files remain in the Zenodo release archive to keep this GitHub repository lighter.

Figure source-data tables are included in `tables/Figure1_source_data.tsv` through `tables/Figure6_source_data.tsv`.

## Example Commands

Rename GFF3/GTF/FASTA records with the v1 mapping:

```bash
python3 scripts/rename_ttenax_annotation_ids.py \
  --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
  --input old_annotation.gff3 \
  --output annotation/Ttenax_annotation_v1.gff3
```

Apply the mapping to a TSV table:

```bash
python3 scripts/apply_gene_id_mapping_to_tables.py \
  --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
  --input old_table.tsv \
  --output mapped_table.tsv \
  --add-old-gene-id
```

Rebuild the multi-omics matrix:

```bash
python3 scripts/build_multiomics_gene_evidence_matrix.py \
  --mapping annotation/Ttenax_gene_id_mapping_v1.tsv \
  --drs-counts DRS_counts.s0.txt \
  --sr-cdna-counts SR_cDNA_counts.s0.txt \
  --lr-cdna-counts LR_cDNA_counts.s0.txt \
  --proteomics-gene-level proteomics_gene_level.tsv \
  --output annotation/Ttenax_multiomics_gene_evidence_matrix_v1.tsv
```

Render the Figure 1 workflow schematic from the plotting scripts:

```bash
Rscript scripts/figure_plotting_scripts/plot_figure1_workflow.R results/Tt_Figure1_workflow_v7
```

## License

The included release data license is CC-BY-4.0. See `LICENSE`.
