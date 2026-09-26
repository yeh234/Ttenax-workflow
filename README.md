# Ttenax-workflow

Compact workflow and release-support repository for the *Trichomonas tenax* genome annotation and Scientific Data revision.

The current `main` branch includes the September 2026 revision: Figure 1–6 plotting scripts `20260916.1`, their exported source metrics, and Supplementary Tables S1–S7. The historical `v1.0.0` tag and [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.19811482) are unchanged. This revision snapshot is not a finalized Zenodo v1.1 release.

This repository intentionally keeps the v1 annotation release files, command-level SOPs, and small scripts and tables needed to document the annotation naming, multi-omics evidence matrix, and manuscript figure source workflows. Large raw or intermediate data such as POD5/FASTQ/BAM, STAR indexes, Flye directories, BUSCO directories, GALBA output, EggNOG output, and InterProScan output are not stored here.

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
  Supplementary_Tables_S1-S7_20260926.xlsx
  Figure1_source_data.tsv
  ...
  Figure6_source_data.tsv
  revision_20260916/

docs/
  figure_legend_additions_v20260916.1.md
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

The original v1.0.0 figure source-data tables remain in `tables/Figure1_source_data.tsv` through `tables/Figure6_source_data.tsv`. The metric exports accompanying the revised plotting scripts are in `tables/revision_20260916/`; the versioned filenames distinguish them from the original exports. Figure 1 is a schematic and does not emit a new metric table.

`tables/Supplementary_Tables_S1-S7_20260926.xlsx` replaces the S1–S5 workbook in the current branch. It contains S1–S5, S6 contig status, and S7 file inventory, plus a README sheet. It is based on the author's `Supplementary_Tables_S1-S7_20260916.xlsx` revision snapshot, which remains unchanged. Five S7 data cells and the status/source note were updated to record the 23 September GenBank/WGS release, retain the pending GCA/crosswalk status, and link the complete protein search database to PRIDE PXD084213. Affected row heights were increased to show the updated text. All other workbook values and formulas remain unchanged.

## Repository records

- NCBI BioProject: PRJNA1457677; BioSample: SAMN57469269; SRA Study: SRP694902.
- GenBank/WGS: JCCIQK000000000; version described: JCCIQK010000000. NCBI notified the author of release on 23 September 2026. The corresponding Assembly (GCA) identifier remains pending confirmation for the final revision.
- Proteomics: [PRIDE PXD084213](https://www.ebi.ac.uk/pride/archive/projects/PXD084213), DOI [10.6019/PXD084213](https://doi.org/10.6019/PXD084213). The author has confirmed public release.
- Annotation archive: [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.19811482). A new Zenodo version will be prepared separately; no new Zenodo DOI is claimed here.

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
Rscript /path/to/Ttenax-workflow/scripts/figure_plotting_scripts/plot_figure1_workflow.R \
  results/Figure1_20260916_v1
```

Run Figure 2–6 from the full analysis project root, not from a bare repository checkout. Their inputs are analysis files under paths such as `flye_len15k/`, `dorado_polish_gpu1/`, `Express/`, and `final_Tt/`, as specified in each script. Those large analysis directories are intentionally excluded from GitHub; the exported metrics are provided for inspection. All six scripts accept an output-prefix argument and use base R with `grid`. The plotting-script version is `20260916.1`.

The revision removes the large figure titles and top summaries, relocates interpretive text into the manuscript legend companion, and reduces unused panel margins. These presentation changes do not change the source calculations.

## License

The included release data license is CC-BY-4.0. See `LICENSE`.
