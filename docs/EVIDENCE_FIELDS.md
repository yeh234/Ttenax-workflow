# Gene-level evidence fields

## Release v1.1 / public GenBank sequence authority

The current primary annotation and sequences follow GCA_061597355.1 / JCCIQK010000000. `GenBank_Product_Name`, `Length_aa`, public accession fields and `Sequence_version` are the current authoritative records. `Original_Analysis_Length_aa`, `Original_Reference_GFF3_Product_Name`, `Original_name_agreement` and `Functional_Product_Name` describe the historical analysis. They are not newly computed public-version annotations.

All 1,270 historical-to-public protein changes are explicitly classified in `Sequence_change_from_analysis`. The original eight partial-CDS reading-frame cases now use the public proteins. Existing functional predictions and RNA/proteomics evidence were not rerun against these proteins. `Functional_annotation_provenance` and `Evidence_provenance` retain this distinction for every row. Proteomics support is historical gene-level evidence, not a assertion that each peptide was revalidated against the current sequence.

STAR model summaries were rebuilt for public exon coordinates while retaining exact coordinate-matched historical observations. The public annotation has 8,446 model introns; the 54 omitted historical boundary introns had zero read support. This is not a new alignment analysis.

## Historical gene-level measurements

The comprehensive table joins the unchanged v1 gene-evidence matrix by AC16WH GeneID. The same gene-level values are repeated for each isoform. They do not establish isoform-specific support.

`Gene_DRS_count`, `Gene_SR_cDNA_count`, and `Gene_LR_cDNA_count` retain the assigned counts from the respective featureCounts inputs. Count units follow their assay and counting protocol; counts are not cross-assay normalized expression values. `*_support` means at least one assigned count; `*_ge5` and `*_ge10` preserve the original count thresholds. RNA_platforms_supported counts supported RNA platforms (0–3); Any_RNA_support and All_three_RNA_support describe detection across those platforms.

Proteomics_support denotes detection in the deposited single-run experiment. Max_unique_peptides and the other proteomics quantities retain the gene-level aggregations in the source matrix; sums of reported peptide counts are not necessarily distinct peptide unions. High_confidence_proteomics_2_unique_peptides, Strong_proteomics_3_unique_peptides and Very_strong_proteomics_5_unique_peptides retain the original thresholds of at least 2, 3 and 5 unique peptides, respectively. Best_protein_by_unique_peptides is a source-matrix selection, not evidence that the other isoforms were individually detected.

Evidence_class preserves the existing category and its priority order: all three RNA platforms plus proteomics with at least two unique peptides; all three RNA plus proteomics; proteomics-supported models; all three RNA; two or more RNA; one RNA; no detected omics support. The legacy label `RNA_and_proteomics_supported` is assigned by the original code's proteomics branch; read the actual RNA_platforms_supported field rather than inferring RNA support from the label alone. Multiomics_support retains the source flag.

No detected evidence means no support in these assays under the tested conditions. It does not mean biological absence, a failed gene prediction or a confirmed lack of expression in other conditions.
