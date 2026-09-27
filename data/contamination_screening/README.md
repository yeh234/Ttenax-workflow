# Trichomonas tenax Hs-4:NIH contamination screening

## Input assembly

- File: `Ttenax_Hs4_NIH_ATCC30207_AC16WH_v1.fsa`
- Assembly state: final two-round annotation reference
- Contigs: 178
- Total length: 99,194,240 bp
- SHA-256: `3eb875ffbb217251b40dc10deb5e1e5ff52cc2bbb37713b23f13aba6092b3062`
- Galaxy history: `Ttenax_Hs4_NCBI_FCS_screening_20260821`

## Screening runs

### NCBI FCS-GX

- Galaxy tool version: 0.5.5+galaxy2
- FCS-GX revision reported in output: `v0.5.5-9-gfc836c1`
- Database: Complete GX database; build date 2023-01-24; 3,025,824 sequences; 709.264 Gbp
- Taxonomic identifier: 43075 (`Trichomonas tenax`)
- Asserted, inferred, and corrected primary division: `prst:monads`
- Aggregate coverage: 0.726714
- Taxonomy records: 180 split-sequence records spanning 99,194,240 bp
- Result classes: 127 `primary-div`, 51 `repeat`, 1 `low-coverage`, and 1 `inconclusive`
- Assigned final divisions: 174 `prst:monads`, 5 `none`, and 1 `anml:insects`
- The only non-primary inconclusive assignment was a 2,055-bp internal segment of `contig_197` (`contig_197~~143650..145704`) with `anml:insects` division coverage of 50%. It was flanked by `prst:monads` segments of the same contig.
- The FCS-GX action report contains zero action rows; no `EXCLUDE`, `TRIM`, or other corrective action was recommended.

### NCBI FCS-adaptor

- Tool version: 0.5.0
- Mode: Eukaryotes
- The report contains only its header and zero adaptor/vector calls.
- The cleaned output retained 178 sequences.

## Interpretation

Neither screen recommended removal or trimming of any contig. The submitted assembly was therefore not altered on the basis of these screens. The small insect-like interval is reported transparently as an inconclusive assignment, not as an actionable contaminant.

## Archived files and checksums

- `FCS_GX_taxonomy_report.tsv`: `28c04e4de0c2d7730793e458c60a166360f54a24eb9a5f9ecbd17c2c0c7ee99b`
- `FCS_GX_action_report.tsv`: `7134aa1a999f07ede1e2f799fdebedcdf7498d2c6d8241ea8d92b21e6cc7c271`
- `FCS_adaptor_report.tsv`: `ef1e9450d34b0d3a110654000dfff742cb3ba9fc697c2ac8c788e0baba3284dd`

