# SOP: Assembly and annotation workflow for the Trichomonas tenax genome project

Project: Trichomonas tenax genome assembly, polishing, and annotation

Repository note: released annotation files are under `annotation/`, workflow scripts are under `scripts/`, and supplementary tables are under `tables/`. Large raw and intermediate files referenced below are external inputs and are not stored in this repository.

Companion figure script:

- `scripts/figure_plotting_scripts/plot_figure1_workflow.R`

Purpose: Provide the final command-level workflow corresponding to Figure 1, from PromethION long-read data through assembly comparison, depth diagnosis, haplotig reduction, polishing, completeness assessment, comparative homology mapping, structural annotation, functional annotation, final annotation merging, and figure rendering.

---

## 0. Core concept of the final workflow

Figure 1 separates three reference states that should not be conflated:

```text
1. Structural backbone assembly
   flye_len15k/assembly.fasta
   101.43 Mb, 212 contigs, N50 924.5 kb

2. Archived single-round polished reference
   final_ref/Ttenax.nohap.dorado_polished.fasta.gz
   99.19 Mb, 178 contigs, N50 939.5 kb

3. Final local annotation-oriented working reference
   dorado_polish_gpu1/draft.polish2.bs6.fasta
   99.19 Mb, 178 contigs
```

The structural backbone was retained for repeat-collapse and exploratory genome-structure analyses. The two-round polished reference was used for comparative homology mapping, GALBA annotation, gffread extraction, EggNOG, T. vaginalis rescue, InterProScan, and final annotation integration.

---

## 1. Directory setup

Use the project workspace that contains raw reads, intermediate assemblies, and annotation outputs.

Define portable locations for the project, external references, containers, and
raw POD5 input, then run the workflow from the project workspace:

```bash
PROJECT="${PROJECT:-/path/to/Tt_Genome}"
REFERENCE_ROOT="${REFERENCE_ROOT:-/path/to/reference}"
CONTAINER_DIR="${CONTAINER_DIR:-/path/to/containers}"
POD5_DIR="${POD5_DIR:-/path/to/Pod5}"

cd "$PROJECT"
```

Suggested directory layout:

```bash
mkdir -p \
  Raw \
  dorado_sup_fastq \
  dorado_sup_moves \
  flye_hq flye_raw \
  flye_asm60x_hq flye_asm60x_raw \
  flye_len15k flye_len20k \
  flye_longbiased_8p22Gb_hq \
  final_ref \
  dorado_polish_gpu1 \
  busco \
  galba_out \
  eggnogHit \
  interproscan_results \
  final_Tt \
  figures
```

Useful variable block:

```bash
THREADS=84
RAW_FASTQ="$PROJECT/Raw/Tt_GenomeSUP.fastq"
GENOME_SIZE=140m
TV_PROT="$REFERENCE_ROOT/TrichDBv68/TrichDB-68_TvaginalisG32022_AnnotatedProteins.fasta"
```

---

## 2. Dorado SUP basecalling

Input:

```text
$POD5_DIR
```

Model:

```text
dna_r10.4.1_e8.2_400bps_sup@v5.2.0
```

### 2.1 Export SUP FASTQ

```bash
./dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 "$POD5_DIR" \
  -r \
  --emit-fastq \
  --emit-summary \
  -o dorado_sup_fastq
```

### 2.2 Export move-aware BAM for polishing

```bash
./dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 "$POD5_DIR" \
  -r \
  --emit-summary \
  --emit-moves \
  --output-dir dorado_sup_moves
```

### 2.3 Concatenate FASTQ if needed

```bash
cat dorado_sup_fastq/*.fastq > Raw/Tt_GenomeSUP.fastq
```

### 2.4 Sequencing yield check

```bash
seqkit stats Raw/Tt_GenomeSUP.fastq
```

Final dataset used in Figure 1:

```text
Reads: 3,726,422
Bases: 11,161,804,683 bp
```

---

## 3. Flye assembly comparison

The goal was to compare several read-subset and read-mode strategies, then select the best structural backbone.

### 3.1 Full-read assemblies

```bash
flye --nano-hq Raw/Tt_GenomeSUP.fastq \
  -g 140m \
  -o flye_hq \
  -t 84

flye --nano-raw Raw/Tt_GenomeSUP.fastq \
  -g 140m \
  -o flye_raw \
  -t 84
```

Observed summary:

```text
flye_hq:  77.6 Mb, N50 245 kb
flye_raw: 62.0 Mb, N50 170 kb
```

### 3.2 asm60x assemblies

```bash
flye --nano-hq asm60x.filtlong.fastq.gz \
  -g 140m \
  -o flye_asm60x_hq \
  -t 84

flye --nano-raw asm60x.filtlong.fastq.gz \
  -g 140m \
  -o flye_asm60x_raw \
  -t 84
```

Observed summary:

```text
flye_asm60x_hq:  96.9 Mb, N50 734 kb, 331 contigs
flye_asm60x_raw: 67.9 Mb, N50 214 kb, 834 contigs
```

### 3.3 Length-filtered assemblies

Extract reads >=15 kb.

```bash
seqkit seq -m 15000 Raw/Tt_GenomeSUP.fastq | \
  gzip -c > Tt_hq.len15k.fastq.gz
```

Assemble reads >=15 kb.

```bash
flye --nano-hq Tt_hq.len15k.fastq.gz \
  -g 140m \
  -o flye_len15k \
  -t 84
```

Extract reads >=20 kb.

```bash
seqkit seq -m 20000 Raw/Tt_GenomeSUP.fastq | \
  gzip -c > Tt_hq.len20k.fastq.gz
```

Assemble reads >=20 kb.

```bash
flye --nano-hq Tt_hq.len20k.fastq.gz \
  -g 140m \
  -o flye_len20k \
  -t 84
```

Observed summary:

```text
flye_len15k: 101.43 Mb, 212 contigs, N50 924.9 kb, mean coverage ~31x
flye_len20k: 102.06 Mb, 273 contigs, N50 754 kb, mean coverage ~19x
```

### 3.4 Long-biased assembly

Build concatenated long-biased input.

```bash
cat reads.ge15k.fastq.gz \
  reads.10to15k.fastq.gz \
  reads.8to10k.fastq.gz \
  reads.6to8k.fastq.gz \
  reads.4to6k.fastq.gz \
  reads.3to4k.fastq.gz \
  > asm_longbiased_8p22Gb.fastq.gz
```

Run Flye.

```bash
flye --nano-hq asm_longbiased_8p22Gb.fastq.gz \
  -g 140m \
  -o flye_longbiased_8p22Gb_hq \
  -t 84
```

Observed summary:

```text
flye_longbiased_8p22Gb_hq: 96.7 Mb, 320 contigs, N50 ~733 kb
```

### 3.5 Backbone decision

Retain:

```text
flye_len15k/assembly.fasta
```

Reason:

```text
101.43 Mb, 212 contigs, N50 924.5-924.9 kb
Better contiguity and effective support than the len20k and long-biased assemblies.
```

---

## 4. Read remapping and depth-based diagnosis

### 4.1 Remap all reads to the selected backbone

```bash
minimap2 -t 84 -ax map-ont flye_len15k/assembly.fasta Raw/Tt_GenomeSUP.fastq | \
  samtools sort -@ 16 -o tt.len15k.allreads.map.bam

samtools index tt.len15k.allreads.map.bam
```

### 4.2 Per-contig mean depth

```bash
samtools depth -a tt.len15k.allreads.map.bam | \
  awk '{cov[$1]+=$3; len[$1]++} END{for (c in cov) printf "%s\t%.2f\t%d\n", c, cov[c]/len[c], len[c]}' | \
  sort -k2,2n > contig_mean_depth.tsv
```

Observed depth statistics:

```text
Representative single-copy depth proxy: 102.56x
90th percentile depth: 119.01x
```

### 4.3 Genome size estimate

Formula:

```text
estimated genome size = total basecalled bases / representative single-copy depth
```

Calculation:

```bash
awk 'BEGIN{printf "%.0f\n", 11161804683 / 102.56}'
```

Observed estimate:

```text
108,831,949 bp
108.83 Mb
```

### 4.4 Depth-based contig classes

Representative local categories:

```text
very-low-depth contigs: 3
haplotig-like contigs: 34
collapsed-repeat contigs using >2C: 6
repeat_gt1p5C contigs using >1.5C: 11
```

Example command pattern for threshold-based extraction from a depth table:

```bash
awk '$2 < 20 {print $0}' contig_mean_depth.tsv > verylow.tsv
awk '$2 >= 20 && $2 < 70 {print $0}' contig_mean_depth.tsv > haplotig_like.tsv
awk '$2 > 205.12 {print $0}' contig_mean_depth.tsv > collapsed.repeat.tsv
awk '$2 > 153.84 {print $0}' contig_mean_depth.tsv > repeat_gt1p5C.tsv
```

Note: thresholds above reflect the local 102.56x single-copy estimate and should be regenerated if the depth estimate changes.

---

## 5. Haplotig reduction

### 5.1 Extract haplotig-like contig IDs

```bash
cut -f1 haplotig_like.tsv > haplotig.ids
```

### 5.2 Build the haplotig-reduced draft

```bash
seqkit grep -v -f haplotig.ids \
  flye_len15k/assembly.fasta \
  > flye_len15k/assembly.no_haplotig.fasta
```

Observed product:

```text
flye_len15k/assembly.no_haplotig.fasta
178 contigs
99.22 Mb
```

This haplotig-reduced assembly became the polishing input.

---

## 6. Single-round Dorado polishing and archived reference

This produced the archived single-round reference state.

### 6.1 Align move-aware BAMs to the no-haplotig draft

```bash
dorado aligner assembly.no_haplotig.fasta PBE76679.bam > PBE76679.aln.bam
dorado aligner assembly.no_haplotig.fasta PBG38932.bam > PBG38932.aln.bam
```

### 6.2 Merge, sort, and index

```bash
samtools merge -@ 16 -o aln.merged.bam PBE76679.aln.bam PBG38932.aln.bam
samtools sort -@ 16 -o aln.merged.sorted.bam aln.merged.bam
samtools index aln.merged.sorted.bam
```

### 6.3 Polish once

```bash
dorado polish -t 40 --infer-threads 20 -x cpu \
  --ignore-read-groups \
  aln.merged.sorted.bam assembly.no_haplotig.fasta \
  > Tt.nohap.polished.fasta
```

Observed product:

```text
Tt.nohap.polished.fasta
178 contigs
99,190,277 bp
N50 939,460 bp
Longest contig 3,767,371 bp
```

### 6.4 Archive the single-round reference

```bash
mkdir -p final_ref
cp Tt.nohap.polished.fasta final_ref/Ttenax.nohap.dorado_polished.fasta
gzip -c final_ref/Ttenax.nohap.dorado_polished.fasta > final_ref/Ttenax.nohap.dorado_polished.fasta.gz
```

---

## 7. Two-round Dorado accuracy polishing

This later polishing stage started from the archived single-round reference and generated the final local working reference.

### 7.1 Setup

```bash
SIF="$CONTAINER_DIR/dorado.sif"
MODELDIR="$PROJECT/dorado_models"
READS="$PROJECT/Raw/Tt_GenomeSUP.fastq"
DRAFT0="$PROJECT/final_ref/Ttenax.nohap.dorado_polished.fasta"
WORK="$PROJECT/dorado_polish_gpu1"
TALIGN=48
TPOLISH=24
INFER=2
BS=6

mkdir -p "$WORK"
cd "$WORK"
```

### 7.2 Round 1

```bash
singularity exec --nv "$SIF" dorado aligner \
  -t "$TALIGN" \
  "$DRAFT0" "$READS" | \
  samtools sort -@ "$TPOLISH" -o aln.r1.dorado.bam -

samtools index -@ "$TPOLISH" aln.r1.dorado.bam

samtools view -@ "$TPOLISH" -b -F 0x900 aln.r1.dorado.bam > aln.r1.dorado.primary.bam
samtools index -@ "$TPOLISH" aln.r1.dorado.primary.bam

singularity exec --nv "$SIF" dorado polish \
  -t "$TPOLISH" \
  --infer-threads "$INFER" \
  --batchsize "$BS" \
  --models-directory "$MODELDIR" \
  --ignore-read-groups \
  aln.r1.dorado.primary.bam "$DRAFT0" \
  > draft.polish1.bs${BS}.fasta
```

### 7.3 Round 2

```bash
singularity exec --nv "$SIF" dorado aligner \
  -t "$TALIGN" \
  draft.polish1.bs${BS}.fasta "$READS" | \
  samtools sort -@ "$TPOLISH" -o aln.r2.dorado.bam -

samtools index -@ "$TPOLISH" aln.r2.dorado.bam

samtools view -@ "$TPOLISH" -b -F 0x900 aln.r2.dorado.bam > aln.r2.dorado.primary.bam
samtools index -@ "$TPOLISH" aln.r2.dorado.primary.bam

singularity exec --nv "$SIF" dorado polish \
  -t "$TPOLISH" \
  --infer-threads "$INFER" \
  --batchsize "$BS" \
  --models-directory "$MODELDIR" \
  --ignore-read-groups \
  aln.r2.dorado.primary.bam draft.polish1.bs${BS}.fasta \
  > draft.polish2.bs${BS}.fasta
```

Observed product:

```text
draft.polish1.bs6.fasta: 178 contigs, 99,198,379 bp
draft.polish2.bs6.fasta: 178 contigs, 99,194,240 bp
```

---

## 8. Polishing evaluation

### 8.1 Mapping summaries

```bash
samtools flagstat aln.r1.dorado.primary.bam > aln.r1.dorado.primary.flagstat.txt
samtools flagstat aln.r2.dorado.primary.bam > aln.r2.dorado.primary.flagstat.txt

cat aln.r1.dorado.primary.flagstat.txt
cat aln.r2.dorado.primary.flagstat.txt
```

Observed:

```text
Round 1 primary reads mapped: 3,666,839 / 3,726,422 = 98.40%
Round 2 primary reads mapped: 3,666,845 / 3,726,422 = 98.40%
```

### 8.2 Residual variant calling pattern

Command pattern:

```bash
bcftools mpileup -Ou -f draft.polish1.bs6.fasta aln.r1.dorado.primary.bam | \
  bcftools call -mv -Oz -o polish1.raw.vcf.gz

bcftools mpileup -Ou -f draft.polish2.bs6.fasta aln.r2.dorado.primary.bam | \
  bcftools call -mv -Oz -o polish2.raw.vcf.gz

bcftools index polish1.raw.vcf.gz
bcftools index polish2.raw.vcf.gz
```

Strict ORF-oriented filter:

```bash
bcftools view -i 'QUAL>=30 && INFO/DP>=20' polish1.raw.vcf.gz > polish1.strict.Q30.DP20.vcf
bcftools view -i 'QUAL>=30 && INFO/DP>=20' polish2.raw.vcf.gz > polish2.strict.Q30.DP20.vcf
```

Observed strict comparison:

```text
Strict HQ INDEL count: polish1 2,358 -> polish2 2,306
Strict HQ INDEL total bp: polish1 31,295 -> polish2 26,323
```

Decision:

```text
Use draft.polish2.bs6.fasta as the final local annotation-oriented working reference.
```

---

## 9. BUSCO and MetaEuk completeness assessment

### 9.1 Default BUSCO

```bash
busco -i draft.polish2.bs6.fasta \
  -m genome \
  -l eukaryota_odb12 \
  -c 48 \
  -o busco_odb12_p2
```

Observed local summary:

```text
C: 31.8%
S: 31.0%
D: 0.8%
F: 2.3%
M: 65.9%
n: 129
```

### 9.2 MetaEuk BUSCO

```bash
busco -i draft.polish1.bs6.fasta \
  -m genome \
  -l eukaryota_odb12 \
  -c 48 \
  -o busco_odb12_p1_metaeuk \
  --metaeuk

busco -i draft.polish2.bs6.fasta \
  -m genome \
  -l eukaryota_odb12 \
  -c 48 \
  -o busco_odb12_p2_metaeuk \
  --metaeuk

busco -i draft.polish2.bs6.fasta \
  -m genome \
  --auto-lineage-euk \
  --metaeuk \
  -c 48 \
  -o busco_p2_auto_euk_metaeuk \
  -f
```

Observed MetaEuk summary:

```text
C: 40.3% (52/129)
S: 37.2% (48/129)
D: 3.1% (4/129)
F: 13.2% (17/129)
M: 46.5% (60/129)
```

Interpretation: BUSCO was treated as a supplementary completeness indicator for this lineage.

---

## 10. Comparative homology mapping with miniprot

### 10.1 Align T. vaginalis proteins to the polished T. tenax genome

```bash
miniprot -t 48 --gff draft.polish2.bs6.fasta \
  "$TV_PROT" \
  > Tt_alignmentTrichDBv68.gff
```

### 10.2 Summarize mapped Tv protein targets

```bash
awk '$3=="mRNA"' Tt_alignmentTrichDBv68.gff | wc -l

awk '$3=="mRNA"{
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^Target=/){
      sub(/^Target=/,"",a[i]);
      split(a[i],b," ");
      print b[1]
    }
  }
}' Tt_alignmentTrichDBv68.gff | sort -u | wc -l
```

Observed:

```text
Total mRNA alignment records: 127,253
Unique Tv protein targets with at least one genomic hit: 42,983
Tv reference proteome denominator: 72,290
Mapped Tv target fraction: 59.5%
```

---

## 11. GALBA structural gene annotation

### 11.1 Setup

```bash
GENOME="$PROJECT/dorado_polish_gpu1/draft.polish2.bs6.fasta"
PROTEIN="$TV_PROT"
GALBA_SIF="$CONTAINER_DIR/galba.sif"
GMES_HOST=/path/to/gmes
```

### 11.2 Run GALBA

```bash
singularity exec \
  -B ${GMES_HOST}:/opt/gmes \
  -B "$PROJECT" \
  -B "$REFERENCE_ROOT" \
  ${GALBA_SIF} \
  bash -c 'export PATH=/opt/gmes:$PATH; galba.pl --genome='"${GENOME}"' --prot_seq='"${PROTEIN}"' --threads=84 --crf --gff3 --workingdir=galba_out'
```

### 11.3 Count structural features

```bash
grep -P "\tgene\t" galba_out/galba.gtf | wc -l
grep -P "\tmRNA\t" galba_out/galba.gff3 | wc -l
```

Observed:

```text
GALBA gene features: 33,725
GALBA mRNA features: 33,845
```

---

## 12. Extract predicted CDS and protein sequences

### 12.1 Run gffread

```bash
gffread galba_out/galba.gff3 \
  -g draft.polish2.bs6.fasta \
  -y Tt_predicted_proteins.aa.fasta \
  -w Tt_predicted_cds.nt.fasta
```

### 12.2 Count extracted sequences

```bash
grep ">" -c Tt_predicted_proteins.aa.fasta
grep ">" -c Tt_predicted_cds.nt.fasta
```

Observed:

```text
Predicted proteins: 33,845
Predicted CDS sequences: 33,845
```

---

## 13. EggNOG functional annotation

### 13.1 Clean the predicted protein FASTA

```bash
cp Tt_predicted_proteins.aa.fasta Tt_predicted_proteins.aa.fasta.bak
sed -i 's/\.//g' Tt_predicted_proteins.aa.fasta
sed -i 's/\*//g' Tt_predicted_proteins.aa.fasta
```

### 13.2 Run EggNOG-mapper

```bash
emapper.py -i Tt_predicted_proteins.aa.fasta \
  --data_dir /dev/shm/eggnog_data \
  -o Tt_eggnog_annotation_final \
  --cpu 48 \
  -m diamond \
  --no_file_comments \
  --override
```

Observed:

```text
Total proteins: 33,845
EggNOG-annotated protein IDs: 19,338
EggNOG-unannotated protein IDs: 14,507
Initial EggNOG annotation rate: 57.1%
```

---

## 14. T. vaginalis BLAST transfer and EggNOG-unannotated rescue

### 14.1 Build the Tv protein BLAST database

```bash
mkdir -p Tv_prot_db

makeblastdb \
  -in "$TV_PROT" \
  -dbtype prot \
  -out Tv_prot_db/Tv_prot_db
```

### 14.2 BLAST the full T. tenax protein set to Tv

```bash
blastp \
  -query Tt_predicted_proteins.aa.fasta \
  -db Tv_prot_db/Tv_prot_db \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
  -evalue 1e-10 \
  -max_target_seqs 1 \
  -num_threads 48 \
  -out Tt_vs_Tv_annotation.txt
```

Observed:

```text
Full Tv transfer rows: 34,849
Unique Tt query IDs represented: 33,594
```

### 14.3 Rescue the EggNOG-unannotated subset

Prepare ID lists.

```bash
grep "^>" Tt_predicted_proteins.aa.fasta | sed 's/^>//' | cut -d' ' -f1 | sort > all_ids.txt

awk 'BEGIN{FS="\t"} $1 !~ /^#/ {print $1}' \
  eggnogHit/Tt_eggnog_annotation_final.emapper.annotations | \
  sort -u > annotated_ids.txt

comm -23 all_ids.txt annotated_ids.txt > unannotated_ids.txt
```

Extract unannotated proteins.

```bash
seqkit grep -f unannotated_ids.txt \
  Tt_predicted_proteins.aa.fasta \
  > unannotated.fasta
```

Run rescue BLAST.

```bash
blastp \
  -query unannotated.fasta \
  -db Tv_prot_db/Tv_prot_db \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads 48 \
  -max_target_seqs 1 \
  -out unannotated_vs_Tv.txt
```

Add Tv descriptions.

```bash
python add_description.py
```

Observed:

```text
EggNOG-unannotated input proteins: 14,507
Rescue-stage rows in unannotated_vs_Tv_with_desc.txt: 15,047
Unique rescue query IDs: 14,459
```

---

## 15. InterProScan domain validation

### 15.1 Run InterProScan

Command pattern:

```bash
interproscan.sh \
  -i Tt_predicted_proteins.aa.fasta \
  -f TSV,GFF3,XML \
  -goterms \
  -iprlookup \
  -pa \
  -cpu 48 \
  -o interproscan_results/Tt_predicted_proteins.aa.fasta.tsv
```

If using a container, bind the working directory and InterProScan data/cache directories as needed.

### 15.2 Parse InterProScan output

```bash
python parse_ips.py
```

Observed parsed summary:

```text
Unique proteins with at least one InterProScan record: 27,223
BspA-containing proteins: 319
Additional LRR-containing proteins: 759
Armadillo-containing proteins: 22
Kinase-containing proteins: 1,293
```

### 15.3 Resolved status of the original EggNOG-unannotated fraction

Observed:

```text
EggNOG-unannotated proteins: 14,507
With Tv rescue hit: 14,459
With InterProScan record: 10,156
With Tv rescue and/or InterProScan support: 14,484
Without Tv rescue and without InterProScan support: 23
```

---

## 16. Merge structural and functional annotation evidence

### 16.1 Merge GALBA, EggNOG, Tv rescue, and InterProScan into final GFF3

The local merge workflow used these input roles:

```python
base_gff = "galba_out/galba.gff3"
output_gff = "Tt_Final_Annotation.gff3"
eggnog_file = "eggnogHit/Tt_eggnog_annotation_final.emapper.annotations"
blast_file = "eggnogHit/unannotated_vs_Tv_with_desc.txt"
ips_file = "interproscan_results/Tt_predicted_proteins.aa.fasta.tsv"
```

Run:

```bash
python merge_galba_annotation.py
```

Recorded merge message:

```text
Success! Annotated 33822 mRNAs in Tt_Final_Annotation.gff3
```

Validate final mRNA count:

```bash
grep -P "\tmRNA\t" Tt_Final_Annotation.gff3 | wc -l
```

Expected:

```text
33,845
```

### 16.2 Build the comprehensive annotation table

```bash
python merge_tsv_fix_length.py
```

Validate row count:

```bash
wc -l Tt_Comprehensive_Annotation_Table.tsv
```

Expected:

```text
33,846 rows including header
33,845 transcript rows
```

Observed final naming-source distribution:

```text
Tv_Rescue: 14,201
EggNOG: 13,583
InterProScan: 2,554
None: 3,507
```

Observed functional coverage:

```text
GO terms: 13,251 proteins
InterPro IDs: 19,538 proteins
KEGG pathway assignments: 7,928 proteins
```

---

## 17. Build and verify the released annotation package

### 17.1 Copy final outputs

```bash
mkdir -p final_Tt

cp dorado_polish_gpu1/draft.polish2.bs6.fasta final_Tt/
cp dorado_polish_gpu1/Tt_Final_Annotation.gff3 final_Tt/
cp dorado_polish_gpu1/Tt_predicted_proteins.aa.fasta final_Tt/
cp dorado_polish_gpu1/Tt_predicted_cds.nt.fasta final_Tt/
cp dorado_polish_gpu1/Tt_Comprehensive_Annotation_Table.tsv final_Tt/
```

### 17.2 Verify copied files

```bash
cmp dorado_polish_gpu1/Tt_Final_Annotation.gff3 final_Tt/Tt_Final_Annotation.gff3
cmp dorado_polish_gpu1/Tt_predicted_proteins.aa.fasta final_Tt/Tt_predicted_proteins.aa.fasta
cmp dorado_polish_gpu1/Tt_predicted_cds.nt.fasta final_Tt/Tt_predicted_cds.nt.fasta
cmp dorado_polish_gpu1/Tt_Comprehensive_Annotation_Table.tsv final_Tt/Tt_Comprehensive_Annotation_Table.tsv
```

Primary released outputs:

```text
final_Tt/Tt_Final_Annotation.gff3
final_Tt/Tt_predicted_proteins.aa.fasta
final_Tt/Tt_predicted_cds.nt.fasta
final_Tt/Tt_Comprehensive_Annotation_Table.tsv
```

---

## 18. Optional downstream multi-omics validation

The detailed multi-omics command SOP is maintained separately:

```text
SOP/multiomics_validation_SOP.md
```

Core final products:

```text
annotation/Ttenax_multiomics_gene_evidence_matrix_v1.tsv
tables/Figure6_source_data.tsv
figures/Figure6_multiomics_support.png
```

Redraw Figure 6:

```bash
cd "$PROJECT"

Rscript scripts/figure_plotting_scripts/plot_figure6_multiomics_support.R \
  figures/Figure6_multiomics_support
```

Add STAR-supported intron annotations to GFF3/GTF:

```bash
cd "$PROJECT"

python3 scripts/add_star_intron_support_to_annotations.py \
  --gff3 Express/Tt_Final_Annotation.gff3 \
  --gtf Express/Tt_Final_Annotation.gtf \
  --supported Express/intron_check/annotation_introns_supported_by_STAR.tsv \
  --supported-min5 Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv \
  --out-prefix Express/Tt_Final_Annotation.STAR_intron_support
```

Validate STAR-supported GTF:

```bash
rg -c $'\tintron\t' Express/Tt_Final_Annotation.STAR_intron_support.gtf
rg -c 'star_junction_support "exact"' Express/Tt_Final_Annotation.STAR_intron_support.gtf
rg -c 'star_junction_support "no"' Express/Tt_Final_Annotation.STAR_intron_support.gtf
```

Expected:

```text
44
44
0
```

Detailed STAR-supported intron SOP:

```text
SOP/STAR_intron_support_SOP.md
```

---

## 19. Redraw Figure 1

From the project workspace:

```bash
cd "$PROJECT"

mkdir -p figures

Rscript scripts/figure_plotting_scripts/plot_figure1_workflow.R \
  figures/Figure1_workflow
```

Expected outputs:

```text
figures/Figure1_workflow.png
figures/Figure1_workflow.pdf
```

---

## 20. Figure 1 final wording checklist

Use the following wording/logic consistently in manuscript text, figure legends, and SOPs:

```text
1. The len15k Flye assembly is the structural backbone, not the final annotation reference.
2. The no-haplotig single-round Dorado polish is an archived reference state.
3. draft.polish2.bs6.fasta is the final local annotation-oriented working reference.
4. BUSCO/MetaEuk are supplementary completeness metrics for this lineage.
5. GALBA is the realized structural annotation route.
6. Functional annotation is a three-layer workflow: EggNOG, Tv rescue, and InterProScan.
7. The final annotation package contains 33,845 transcript-level protein/CDS records.
8. The initial 57.1% EggNOG annotation rate should not be presented as final annotation coverage.
9. Tv BLAST rows are not equivalent to unique covered T. tenax genes.
10. STAR-supported introns are supportive transcriptomic evidence, not targeted validation.
```

---

## 21. Workflow outputs

Assembly and reference states:

```text
flye_len15k/assembly.fasta
flye_len15k/assembly.no_haplotig.fasta
final_ref/Ttenax.nohap.dorado_polished.fasta.gz
dorado_polish_gpu1/draft.polish1.bs6.fasta
dorado_polish_gpu1/draft.polish2.bs6.fasta
```

Annotation products:

```text
galba_out/galba.gff3
Tt_predicted_proteins.aa.fasta
Tt_predicted_cds.nt.fasta
Tt_Final_Annotation.gff3
Tt_Comprehensive_Annotation_Table.tsv
final_Tt/Tt_Final_Annotation.gff3
final_Tt/Tt_predicted_proteins.aa.fasta
final_Tt/Tt_predicted_cds.nt.fasta
final_Tt/Tt_Comprehensive_Annotation_Table.tsv
```

Evidence products:

```text
Tt_alignmentTrichDBv68.gff
Tt_vs_Tv_annotation.txt
eggnogHit/Tt_eggnog_annotation_final.emapper.annotations
eggnogHit/unannotated_vs_Tv_with_desc.txt
interproscan_results/Tt_predicted_proteins.aa.fasta.tsv
Express/Tt_Final_Annotation.STAR_intron_support.gff3
Express/Tt_Final_Annotation.STAR_intron_support.gtf
Express/multiomics_gene_evidence_matrix.tsv
```

Figure products:

```text
figures/Figure1_workflow.png
figures/Figure1_workflow.pdf
figures/Figure6_multiomics_support.png
```
