# STAR Intron Support

This directory contains the compact STAR-supported intron evidence table for the
*Trichomonas tenax* annotation v1.0.0 release.

Included:

```text
Ttenax_STAR_supported_introns_v1.tsv
```

The table reports annotated introns with exact STAR junction evidence, including
gene ID, transcript ID, intron ID, contig, coordinates, strand, unique junction
reads, multimapping junction reads, total junction reads, the STAR annotated
junction flag, and a `star_min5_support` flag.

The full STAR-support-enhanced annotation files are intentionally not mirrored
in this compact GitHub workflow repository:

```text
Ttenax_annotation_STAR_intron_support_v1.gff3
Ttenax_annotation_STAR_intron_support_v1.gtf
```

Those full files are available in the Zenodo/release archive:

```text
Ttenax_annotation_release_v1.0.0.tar.gz
```
