# SOP: Multi-omics validation of the *Trichomonas tenax* genome annotation

**Project:** *Trichomonas tenax* genome assembly and annotation validation  
**Purpose:** Reprocess Nanopore direct RNA sequencing (DRS), Nanopore long-read cDNA, Illumina short-read cDNA, and proteomics evidence to support the final *T. tenax* genome annotation.  
**Final annotation:** `Tt_Final_Annotation.gff3` / `Tt_Final_Annotation.gtf`  
**Final genome:** `draft.polish2.bs6.fasta.gz` or uncompressed `draft.polish2.bs6.fasta`  
**Primary counting mode:** unstranded featureCounts, `-s 0`  

Repository note: the compact GitHub layout stores the frozen v1 matrix at `annotation/Ttenax_multiomics_gene_evidence_matrix_v1.tsv` and the rebuild helper at `scripts/build_multiomics_gene_evidence_matrix.py`. Large DRS/LR cDNA/SR cDNA reads and BAMs remain external inputs.

---

## 0. Directory layout

Recommended working directory:

```bash
/NVMe/Tt2026_RNA
```

Suggested subdirectories:

```bash
mkdir -p \
  DRS LRcDNA SR \
  basecall \
  aln qc counts \
  SR_trim SR_aln SR_counts STAR_index \
  LRcDNA_aln LRcDNA_qc LRcDNA_counts \
  intron_check \
  multiomics_matrix/input multiomics_matrix/output
```

---

## 1. Software environment

### 1.1 Conda / mamba environment for RNA-seq processing

```bash
mamba create -n srna -c conda-forge -c bioconda \
  fastp star subread samtools seqkit minimap2 nanoplot python -y

conda activate srna
```

Confirmed software versions used in this project:

```text
fastp 1.3.2
STAR 2.7.11b
featureCounts 2.1.1
samtools 1.23.1
seqkit 2.13.0
Dorado 0.9.6+0949eb8d
```

### 1.2 POD5 conversion environment

The old `nanopore` environment used Python 3.6 and could not install `pod5`. Use a newer Python environment.

```bash
mamba create -n pod5env -c conda-forge -c bioconda python=3.10 pod5 -y
conda activate pod5env

pod5 --version
```

---

## 2. Convert legacy FAST5 to POD5

The original Nanopore folders were separated into `fast5_pass` and `fast5_fail`. For validation and recovery of legacy data, both were converted and merged.

### 2.1 Direct RNA sequencing, RNA002 R9.4.1

```bash
conda activate pod5env

mkdir -p DRS/pod5

pod5 convert fast5 DRS/fast5_pass/*.fast5 \
  --output DRS/Tt_DRS_pass.pod5

pod5 convert fast5 DRS/fast5_fail/*.fast5 \
  --output DRS/Tt_DRS_fail.pod5

pod5 merge DRS/Tt_DRS_pass.pod5 DRS/Tt_DRS_fail.pod5 \
  --output DRS/Tt_DRS_merged.pod5

ls -lh DRS/*.pod5
```

Expected files:

```text
Tt_DRS_pass.pod5
Tt_DRS_fail.pod5
Tt_DRS_merged.pod5
```

### 2.2 Nanopore cDNA, R9.4.1

```bash
conda activate pod5env

mkdir -p LRcDNA/pod5

pod5 convert fast5 LRcDNA/fast5_pass/*.fast5 \
  --output LRcDNA/Tt_cDNA_pass.pod5

pod5 convert fast5 LRcDNA/fast5_fail/*.fast5 \
  --output LRcDNA/Tt_cDNA_fail.pod5

pod5 merge LRcDNA/Tt_cDNA_pass.pod5 LRcDNA/Tt_cDNA_fail.pod5 \
  --output LRcDNA/Tt_cDNA_merged.pod5

ls -lh LRcDNA/*.pod5
```

---

## 3. Dorado basecalling

### 3.1 Check Dorado and GPU

```bash
dorado -v
nvidia-smi
```

In this project:

```text
Dorado: 0.9.6+0949eb8d
GPU: NVIDIA GeForce RTX 2080 Ti, 11 GB VRAM
Driver: 535.247.01
CUDA: 12.2
```

### 3.2 Available R9.4.1 models observed

```text
dna_r9.4.1_e8_fast@v3.4
dna_r9.4.1_e8_hac@v3.3
dna_r9.4.1_e8_sup@v3.3
dna_r9.4.1_e8_sup@v3.6
rna002_70bps_fast@v3
rna002_70bps_hac@v3
```

### 3.3 Direct RNA basecalling

Input: `DRS/Tt_DRS_merged.pod5`  
Model: `rna002_70bps_hac@v3`

```bash
mkdir -p basecall dorado_model logs

# Optional: download model if not already available.
dorado download --model rna002_70bps_hac@v3 --models-directory dorado_model

# Basecall DRS.
dorado basecaller \
  dorado_model/rna002_70bps_hac@v3 \
  DRS/Tt_DRS_merged.pod5 \
  --device cuda:0 \
  --emit-moves \
  > basecall/Tt_DRS_merged.hac.bam \
  2> logs/Tt_DRS_Merged.hac.basecall.log
```

Notes:

- Dorado may warn: `BAM format does not support U, so RNA output files will include T instead of U`.
- Dorado may warn: `Could not determine sequencing Chemistry from read data`; for legacy RNA002 data this can occur and did not prevent basecalling.

Observed DRS basecalling result:

```text
Simplex reads basecalled: 1,469,841
Simplex reads filtered: 11,669
```

### 3.4 Long-read cDNA basecalling

Input: `LRcDNA/Tt_cDNA_merged.pod5`  
Model: `dna_r9.4.1_e8_sup@v3.6`

```bash
mkdir -p basecall dorado_model logs

# Optional: download model if not already available.
dorado download --model dna_r9.4.1_e8_sup@v3.6 --models-directory dorado_model

# Basecall LR cDNA.
dorado basecaller \
  dorado_model/dna_r9.4.1_e8_sup@v3.6 \
  LRcDNA/Tt_cDNA_merged.pod5 \
  --device cuda:0 \
  --emit-moves \
  > basecall/Tt_cDNA_merged.sup.bam \
  2> logs/Tt_cDNA_merged.sup.basecall.log
```

Observed LR cDNA basecalling result:

```text
Simplex reads basecalled: 12,905,233
Simplex reads filtered: 100
```

---

## 4. Convert Dorado BAM to FASTQ and run read QC

### 4.1 Direct RNA FASTQ

```bash
samtools quickcheck basecall/Tt_DRS_merged.hac.bam && echo OK
samtools view -c basecall/Tt_DRS_merged.hac.bam
samtools stats basecall/Tt_DRS_merged.hac.bam > basecall/Tt_DRS_merged.hac.stats.txt

samtools fastq basecall/Tt_DRS_merged.hac.bam | \
  gzip > basecall/Tt_DRS_merged.hac.fastq.gz

seqkit stats basecall/Tt_DRS_merged.hac.fastq.gz

NanoPlot \
  --fastq basecall/Tt_DRS_merged.hac.fastq.gz \
  -o nanoplot_DRS_hac
```

Observed DRS FASTQ QC:

```text
num_seqs: 1,527,015
sum_len: 1,244,336,578 bp
mean length: 814.9 bp
median length: 641 bp
read length N50: 1,064 bp
mean Q: 13.8
median Q: 16.6
>Q10: 93.2%
>Q15: 69.4%
```

### 4.2 Long-read cDNA FASTQ

```bash
samtools quickcheck basecall/Tt_cDNA_merged.sup.bam && echo OK
samtools view -c basecall/Tt_cDNA_merged.sup.bam
samtools stats basecall/Tt_cDNA_merged.sup.bam > basecall/Tt_cDNA_merged.sup.stats.txt

samtools fastq basecall/Tt_cDNA_merged.sup.bam | \
  gzip > basecall/Tt_cDNA_merged.sup.fastq.gz

seqkit stats basecall/Tt_cDNA_merged.sup.fastq.gz

NanoPlot \
  --fastq basecall/Tt_cDNA_merged.sup.fastq.gz \
  -o nanoplot_LRcDNA_sup
```

Observed raw LR cDNA QC:

```text
num_seqs: 12,906,145
sum_len: 10,289,439,958 bp
mean length: 797.3 bp
median length: 681 bp
read length N50: 921 bp
mean Q: 10.8
median Q: 12.8
>Q10: 83.1%
>Q15: 18.2%
max length: 927,617 bp
```

The extremely long, low-quality reads were treated as likely concatemeric/chimeric artifacts and filtered before genome alignment.

---

## 5. Long-read cDNA filtering

Recommended formal LR cDNA analysis dataset:

- Q-score threshold: `Q >= 10`
- minimum length: `200 bp`
- maximum length: `50,000 bp`

```bash
seqkit seq -Q 10 -m 200 -M 50000 \
  basecall/Tt_cDNA_merged.sup.fastq.gz | \
  gzip > basecall/Tt_cDNA_merged.sup.Q10.min200.max50k.fastq.gz

seqkit seq -Q 10 -m 500 -M 50000 \
  basecall/Tt_cDNA_merged.sup.fastq.gz | \
  gzip > basecall/Tt_cDNA_merged.sup.Q10.min500.max50k.fastq.gz

seqkit stats \
  basecall/Tt_cDNA_merged.sup.fastq.gz \
  basecall/Tt_cDNA_merged.sup.Q10.min200.max50k.fastq.gz \
  basecall/Tt_cDNA_merged.sup.Q10.min500.max50k.fastq.gz
```

Observed filtering result:

```text
raw:                  12,906,145 reads; 10.29 Gb; mean 797.3 bp; max 927,617 bp
Q10.min200.max50k:    10,626,042 reads;  8.80 Gb; mean 828.1 bp; max 48,340 bp
Q10.min500.max50k:     8,368,256 reads;  7.98 Gb; mean 954.1 bp; max 48,340 bp
```

---

## 6. Prepare genome FASTA

If the final genome is compressed:

```bash
gunzip -c draft.polish2.bs6.fasta.gz > draft.polish2.bs6.fasta
```

Optional symlink:

```bash
mkdir -p ref
ln -sf ../draft.polish2.bs6.fasta ref/draft.polish2.bs6.fasta
```

---

## 7. Direct RNA alignment and gene assignment

### 7.1 Align DRS to genome

```bash
mkdir -p aln qc counts

minimap2 -t 36 -ax splice -uf -k14 \
  ref/draft.polish2.bs6.fasta \
  basecall/Tt_DRS_merged.hac.fastq.gz | \
samtools sort -@ 12 -m 2G \
  -o aln/Tt_DRS_merged.hac.genome.bam

samtools index aln/Tt_DRS_merged.hac.genome.bam

samtools flagstat aln/Tt_DRS_merged.hac.genome.bam \
  > qc/Tt_DRS_merged.hac.genome.flagstat.txt

samtools stats aln/Tt_DRS_merged.hac.genome.bam \
  > qc/Tt_DRS_merged.hac.genome.stats.txt

cat qc/Tt_DRS_merged.hac.genome.flagstat.txt
```

Observed DRS mapping:

```text
1,640,612 mapped / 1,803,523 total records = 90.97%
```

### 7.2 DRS featureCounts strandedness test

```bash
featureCounts -L --primary -T 16 -s 0 \
  -a Tt_Final_Annotation.gtf \
  -o counts/Tt_DRS_gene_counts.s0.txt \
  aln/Tt_DRS_merged.hac.genome.bam

featureCounts -L --primary -T 16 -s 1 \
  -a Tt_Final_Annotation.gtf \
  -o counts/Tt_DRS_gene_counts.s1.txt \
  aln/Tt_DRS_merged.hac.genome.bam

featureCounts -L --primary -T 16 -s 2 \
  -a Tt_Final_Annotation.gtf \
  -o counts/Tt_DRS_gene_counts.s2.txt \
  aln/Tt_DRS_merged.hac.genome.bam

cat counts/Tt_DRS_gene_counts.s0.txt.summary
cat counts/Tt_DRS_gene_counts.s1.txt.summary
cat counts/Tt_DRS_gene_counts.s2.txt.summary
```

Observed DRS result:

```text
-s 0 Assigned: 1,392,146
-s 1 Assigned: 1,387,784
-s 2 Assigned:    13,251
```

Primary DRS setting: `-s 0`.

### 7.3 DRS supported genes

```bash
awk 'NR>2 && $7>0 {n++} END{print "genes_with_>=1_read:", n}' counts/Tt_DRS_gene_counts.s0.txt
awk 'NR>2 && $7>=5 {n++} END{print "genes_with_>=5_reads:", n}' counts/Tt_DRS_gene_counts.s0.txt
awk 'NR>2 && $7>=10 {n++} END{print "genes_with_>=10_reads:", n}' counts/Tt_DRS_gene_counts.s0.txt

awk 'NR>2 {print $1"\t"$7}' counts/Tt_DRS_gene_counts.s0.txt | \
  sort -k2,2nr | head -20
```

Observed DRS gene support:

```text
>=1 read: 18,574 genes
>=5 reads: 14,302 genes
>=10 reads: 10,740 genes
```

---

## 8. Illumina short-read cDNA processing

Input paired-end FASTQ:

```text
SR/Tt_HVFMLDSX2_L4_R1.fastq.gz
SR/Tt_HVFMLDSX2_L4_R2.fastq.gz
```

### 8.1 Trim and QC with fastp

```bash
mkdir -p SR_trim fastp_reports

fastp \
  -i SR/Tt_HVFMLDSX2_L4_R1.fastq.gz \
  -I SR/Tt_HVFMLDSX2_L4_R2.fastq.gz \
  -o SR_trim/Tt_HVFMLDSX2.trim.R1.fastq.gz \
  -O SR_trim/Tt_HVFMLDSX2.trim.R2.fastq.gz \
  --thread 16 \
  --detect_adapter_for_pe \
  --html fastp_reports/Tt_HVFMLDSX2.fastp.html \
  --json fastp_reports/Tt_HVFMLDSX2.fastp.json

seqkit stats SR_trim/Tt_HVFMLDSX2.trim.R1.fastq.gz SR_trim/Tt_HVFMLDSX2.trim.R2.fastq.gz
```

Observed trimmed reads:

```text
R1: 37,505,899 reads; 5.17 Gb; mean 137.8 bp
R2: 37,505,899 reads; 5.17 Gb; mean 137.8 bp
```

### 8.2 Build STAR genome index

```bash
mkdir -p STAR_index/Tt_genome

STAR \
  --runThreadN 36 \
  --runMode genomeGenerate \
  --genomeDir STAR_index/Tt_genome \
  --genomeFastaFiles ref/draft.polish2.bs6.fasta \
  --sjdbGTFfile Tt_Final_Annotation.gtf \
  --sjdbOverhang 149
```

### 8.3 STAR alignment

Initial `SortedByCoordinate` failed due to STAR BAM sorting temporary file/open-file limits. Robust approach: output unsorted BAM, then sort with samtools.

```bash
mkdir -p SR_aln
ulimit -n 65535
rm -rf SR_aln/Tt_HVFMLDSX2._STARtmp

STAR \
  --runThreadN 20 \
  --genomeDir STAR_index/Tt_genome \
  --readFilesIn SR_trim/Tt_HVFMLDSX2.trim.R1.fastq.gz SR_trim/Tt_HVFMLDSX2.trim.R2.fastq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix SR_aln/Tt_HVFMLDSX2. \
  --quantMode GeneCounts

samtools sort -@ 12 -m 2G \
  -o SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam \
  SR_aln/Tt_HVFMLDSX2.Aligned.out.bam

samtools index SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam

samtools flagstat SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam \
  > SR_aln/Tt_HVFMLDSX2.flagstat.txt

cat SR_aln/Tt_HVFMLDSX2.Log.final.out
cat SR_aln/Tt_HVFMLDSX2.flagstat.txt
```

Observed STAR mapping:

```text
Input reads: 37,505,899
Uniquely mapped reads: 33,525,278
Uniquely mapped reads %: 89.39%
Splices total: 212,581
Annotated splices: 198,061
```

### 8.4 SR cDNA featureCounts strandedness test

```bash
mkdir -p SR_counts

featureCounts -T 16 -p -s 0 \
  -a Tt_Final_Annotation.gtf \
  -o SR_counts/Tt_HVFMLDSX2.s0.txt \
  SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam

featureCounts -T 16 -p -s 1 \
  -a Tt_Final_Annotation.gtf \
  -o SR_counts/Tt_HVFMLDSX2.s1.txt \
  SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam

featureCounts -T 16 -p -s 2 \
  -a Tt_Final_Annotation.gtf \
  -o SR_counts/Tt_HVFMLDSX2.s2.txt \
  SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam

cat SR_counts/Tt_HVFMLDSX2.s0.txt.summary
cat SR_counts/Tt_HVFMLDSX2.s1.txt.summary
cat SR_counts/Tt_HVFMLDSX2.s2.txt.summary
```

Observed SR cDNA assignment:

```text
-s 0 Assigned: 64,030,411
-s 1 Assigned: 63,615,348
-s 2 Assigned:    422,837
```

Primary SR cDNA setting: `-s 0`.

### 8.5 SR cDNA supported genes

```bash
awk 'NR>2 && $7>0 {n++} END{print "genes_with_>=1_fragment:", n}' SR_counts/Tt_HVFMLDSX2.s0.txt
awk 'NR>2 && $7>=5 {n++} END{print "genes_with_>=5_fragments:", n}' SR_counts/Tt_HVFMLDSX2.s0.txt
awk 'NR>2 && $7>=10 {n++} END{print "genes_with_>=10_fragments:", n}' SR_counts/Tt_HVFMLDSX2.s0.txt

awk 'NR>2 {print $1"\t"$7}' SR_counts/Tt_HVFMLDSX2.s0.txt | \
  sort -k2,2nr | head -20
```

Observed SR cDNA gene support:

```text
>=1 fragment: 19,540 genes
>=5 fragments: 18,559 genes
>=10 fragments: 18,253 genes
```

---

## 9. Short-read splice junction and annotated intron cross-check

This analysis was retained as a minor supporting analysis only. No Sanger sequencing or targeted RT-PCR validation was performed.

### 9.1 Extract annotated introns from GFF3

```bash
mkdir -p intron_check

awk -F'\t' '
BEGIN{OFS="\t"}
$3=="intron"{
  tid=""; iid="";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^Parent=/){sub(/^Parent=/,"",a[i]); tid=a[i]}
    if(a[i] ~ /^ID=/){sub(/^ID=/,"",a[i]); iid=a[i]}
  }
  print $1,$4,$5,$7,tid,iid
}' Tt_Final_Annotation.gff3 | sort -k1,1 -k2,2n > intron_check/annotation_introns.tsv

wc -l intron_check/annotation_introns.tsv
head intron_check/annotation_introns.tsv
```

Observed annotated introns:

```text
8,500 introns
```

### 9.2 Convert STAR SJ.out.tab

```bash
awk '
BEGIN{OFS="\t"}
{
  s="."
  if($4==1) s="+"
  else if($4==2) s="-"
  total=$7+$8
  print $1,$2,$3,s,$7,$8,total,$6
}' SR_aln/Tt_HVFMLDSX2.SJ.out.tab | sort -k1,1 -k2,2n > intron_check/star_junctions.tsv

wc -l intron_check/star_junctions.tsv
head intron_check/star_junctions.tsv
```

Observed STAR unique junctions:

```text
4,476 unique junction rows
```

### 9.3 Exact coordinate + strand matching

```bash
awk 'BEGIN{OFS="\t"} {key=$1":"$2":"$3":"$4; print key,$0}' \
  intron_check/annotation_introns.tsv > intron_check/annotation_introns.keyed.tsv

awk 'BEGIN{OFS="\t"} {key=$1":"$2":"$3":"$4; print key,$0}' \
  intron_check/star_junctions.tsv > intron_check/star_junctions.keyed.tsv

join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 intron_check/annotation_introns.keyed.tsv) \
  <(sort -k1,1 intron_check/star_junctions.keyed.tsv) \
  > intron_check/annotation_introns_supported_by_STAR.tsv

join -t $'\t' -v 1 -1 1 -2 1 \
  <(sort -k1,1 intron_check/annotation_introns.keyed.tsv) \
  <(sort -k1,1 intron_check/star_junctions.keyed.tsv) \
  > intron_check/annotation_introns_NOT_supported_by_STAR.tsv

wc -l intron_check/annotation_introns_supported_by_STAR.tsv
wc -l intron_check/annotation_introns_NOT_supported_by_STAR.tsv
```

Observed exact intron support:

```text
44 annotated introns matched short-read splice junctions exactly
42 transcripts had at least one supported intron
34 transcripts had all annotated introns supported
8 transcripts had partial intron support
6,145 intron-containing transcripts had no exact junction support
```

Suggested manuscript wording:

```text
A total of 44 short-read splice junctions matched annotated intron coordinates exactly, corresponding to 42 transcripts; these observations were retained as supporting transcriptomic evidence but were not further emphasized.
```

---

## 10. Long-read cDNA alignment and gene assignment

### 10.1 Align filtered LR cDNA to genome

```bash
mkdir -p LRcDNA_aln LRcDNA_qc LRcDNA_counts

minimap2 -t 36 -ax splice -k14 \
  ref/draft.polish2.bs6.fasta \
  basecall/Tt_cDNA_merged.sup.Q10.min200.max50k.fastq.gz | \
samtools sort -@ 12 -m 2G \
  -o LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

samtools index LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

samtools flagstat LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam \
  > LRcDNA_qc/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.flagstat.txt

samtools stats LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam \
  > LRcDNA_qc/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.stats.txt

cat LRcDNA_qc/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.flagstat.txt
```

Observed LR cDNA mapping:

```text
Filtered input reads: 10,626,042
Primary mapped: 10,185,924 / 10,626,042 = 95.86%
Total mapped: 12,527,277 / 12,967,395 = 96.61%
```

### 10.2 LR cDNA featureCounts strandedness test

```bash
featureCounts -L --primary -T 16 -s 0 \
  -a Tt_Final_Annotation.gtf \
  -o LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt \
  LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

featureCounts -L --primary -T 16 -s 1 \
  -a Tt_Final_Annotation.gtf \
  -o LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s1.txt \
  LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

featureCounts -L --primary -T 16 -s 2 \
  -a Tt_Final_Annotation.gtf \
  -o LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s2.txt \
  LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

cat LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt.summary
cat LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s1.txt.summary
cat LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s2.txt.summary
```

Observed LR cDNA assignment:

```text
-s 0 Assigned: 8,289,823
-s 1 Assigned: 4,089,169
-s 2 Assigned: 4,268,603
```

Primary LR cDNA setting: `-s 0`.

### 10.3 LR cDNA supported genes

```bash
awk 'NR>2 && $7>0 {n++} END{print "LRcDNA_genes_with_>=1_read:", n}' LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt
awk 'NR>2 && $7>=5 {n++} END{print "LRcDNA_genes_with_>=5_reads:", n}' LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt
awk 'NR>2 && $7>=10 {n++} END{print "LRcDNA_genes_with_>=10_reads:", n}' LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt

awk 'NR>2 {print $1"\t"$7}' LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt | \
  sort -k2,2nr | head -20
```

Observed LR cDNA gene support:

```text
>=1 read: 21,143 genes
>=5 reads: 18,145 genes
>=10 reads: 17,277 genes
```

---

## 11. Proteomics gene-level parsing

Primary proteomics report:

```text
20260312-PT-15_100-150min-tt.txt
```

This file contains protein-level `Master Protein` rows and peptide-level rows. Protein IDs use transcript-like IDs such as `g16424t1`, which were converted to gene-level IDs such as `g16424`.

### 11.1 Create parser script

```bash
cat > parse_proteomics_to_gene_level.py <<'PY'
#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def protein_to_gene_id(accession: str) -> str:
    accession = accession.strip()
    accession = re.sub(r't\d+$', '', accession)
    accession = re.sub(r'\.t\d+$', '', accession)
    return accession


def to_float(x):
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "na", "n/a"}:
        return ""
    try:
        return float(x)
    except ValueError:
        return ""


def to_int(x):
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "na", "n/a"}:
        return 0
    try:
        return int(float(x))
    except ValueError:
        return 0


def main():
    parser = argparse.ArgumentParser(
        description="Parse Proteome Discoverer/Mascot-style proteomics txt and convert protein IDs to gene-level IDs."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="proteomics_gene_level.tsv")
    parser.add_argument("--min-unique-peptides", type=int, default=1)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    records = []

    with input_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        required_cols = [
            "Master", "Accession", "Coverage [%]", "# Peptides", "# PSMs",
            "# Unique Peptides", "# AAs", "MW [kDa]", "calc. pI",
            "Score Mascot: Mascot",
        ]
        missing = [c for c in required_cols if c not in reader.fieldnames]
        if missing:
            raise ValueError("Missing required columns: " + ", ".join(missing))

        for row in reader:
            if row.get("Master", "").strip() != "Master Protein":
                continue
            protein_id = row.get("Accession", "").strip()
            if not protein_id:
                continue
            gene_id = protein_to_gene_id(protein_id)
            unique_peptides = to_int(row.get("# Unique Peptides", 0))
            if unique_peptides < args.min_unique_peptides:
                continue
            records.append({
                "GeneID": gene_id,
                "ProteinID": protein_id,
                "Description": row.get("Description", "").strip(),
                "Coverage_percent": row.get("Coverage [%]", "").strip(),
                "Peptides": to_int(row.get("# Peptides", 0)),
                "PSMs": to_int(row.get("# PSMs", 0)),
                "Unique_peptides": unique_peptides,
                "AAs": to_int(row.get("# AAs", 0)),
                "MW_kDa": row.get("MW [kDa]", "").strip(),
                "calc_pI": row.get("calc. pI", "").strip(),
                "Mascot_score": row.get("Score Mascot: Mascot", "").strip(),
                "Protein_groups": row.get("# Protein Groups", "").strip(),
                "Abundance_scaled": row.get("Abundances (Scaled): F337: Sample, n/a", "").strip(),
                "Abundance": row.get("Abundance: F337: Sample, n/a", "").strip(),
                "Found_in_sample": row.get("Found in Sample: [S333] F337: Sample, n/a", "").strip(),
                "Modifications": row.get("Modifications", "").strip(),
            })

    gene_summary = {}

    for r in records:
        gid = r["GeneID"]
        if gid not in gene_summary:
            gene_summary[gid] = {
                "GeneID": gid,
                "ProteinIDs": [],
                "Descriptions": [],
                "Max_coverage_percent": 0.0,
                "Total_peptides_reported": 0,
                "Total_PSMs": 0,
                "Max_unique_peptides": 0,
                "Sum_unique_peptides_reported": 0,
                "Max_Mascot_score": 0.0,
                "Best_protein_by_unique_peptides": "",
                "Proteomics_support": "yes",
                "High_confidence_2_unique_peptides": "no",
                "Strong_support_3_unique_peptides": "no",
                "Very_strong_support_5_unique_peptides": "no",
            }

        g = gene_summary[gid]
        g["ProteinIDs"].append(r["ProteinID"])
        if r["Description"]:
            g["Descriptions"].append(r["Description"])

        cov = to_float(r["Coverage_percent"])
        if cov != "":
            g["Max_coverage_percent"] = max(g["Max_coverage_percent"], cov)

        mascot = to_float(r["Mascot_score"])
        if mascot != "":
            g["Max_Mascot_score"] = max(g["Max_Mascot_score"], mascot)

        g["Total_peptides_reported"] += r["Peptides"]
        g["Total_PSMs"] += r["PSMs"]
        g["Sum_unique_peptides_reported"] += r["Unique_peptides"]

        if r["Unique_peptides"] > g["Max_unique_peptides"]:
            g["Max_unique_peptides"] = r["Unique_peptides"]
            g["Best_protein_by_unique_peptides"] = r["ProteinID"]

    for gid, g in gene_summary.items():
        max_up = g["Max_unique_peptides"]
        if max_up >= 2:
            g["High_confidence_2_unique_peptides"] = "yes"
        if max_up >= 3:
            g["Strong_support_3_unique_peptides"] = "yes"
        if max_up >= 5:
            g["Very_strong_support_5_unique_peptides"] = "yes"
        g["ProteinIDs"] = ";".join(sorted(set(g["ProteinIDs"])))
        g["Descriptions"] = ";".join(sorted(set(g["Descriptions"])))
        g["Max_coverage_percent"] = round(g["Max_coverage_percent"], 3)
        g["Max_Mascot_score"] = round(g["Max_Mascot_score"], 3)

    out_cols = [
        "GeneID", "ProteinIDs", "Descriptions", "Proteomics_support",
        "Max_unique_peptides", "Sum_unique_peptides_reported",
        "Total_peptides_reported", "Total_PSMs", "Max_coverage_percent",
        "Max_Mascot_score", "Best_protein_by_unique_peptides",
        "High_confidence_2_unique_peptides", "Strong_support_3_unique_peptides",
        "Very_strong_support_5_unique_peptides",
    ]

    def gene_sort_key(x):
        m = re.search(r"\d+", x)
        return int(m.group()) if m else 0

    with output_path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_cols)
        writer.writeheader()
        for gid in sorted(gene_summary.keys(), key=gene_sort_key):
            writer.writerow({c: gene_summary[gid][c] for c in out_cols})

    print(f"Input protein-level records retained: {len(records)}")
    print(f"Gene-level records written: {len(gene_summary)}")
    print(f"Output: {output_path}")
    print(f"Genes with >=1 unique peptide: {len(gene_summary)}")
    print(f"Genes with >=2 unique peptides: {sum(1 for g in gene_summary.values() if g['Max_unique_peptides'] >= 2)}")
    print(f"Genes with >=3 unique peptides: {sum(1 for g in gene_summary.values() if g['Max_unique_peptides'] >= 3)}")
    print(f"Genes with >=5 unique peptides: {sum(1 for g in gene_summary.values() if g['Max_unique_peptides'] >= 5)}")


if __name__ == "__main__":
    main()
PY

chmod +x parse_proteomics_to_gene_level.py
```

### 11.2 Run proteomics parser

```bash
python parse_proteomics_to_gene_level.py \
  -i 20260312-PT-15_100-150min-tt.txt \
  -o proteomics_gene_level.tsv

head -5 proteomics_gene_level.tsv
wc -l proteomics_gene_level.tsv
```

### 11.3 Export proteomics-supported gene sets

```bash
awk 'NR>1 && $5>=1 {print $1}' proteomics_gene_level.tsv \
  > proteomics_genes_ge1_unique_peptide.txt

awk 'NR>1 && $5>=2 {print $1}' proteomics_gene_level.tsv \
  > proteomics_genes_ge2_unique_peptides.txt

awk 'NR>1 && $5>=3 {print $1}' proteomics_gene_level.tsv \
  > proteomics_genes_ge3_unique_peptides.txt

awk 'NR>1 && $5>=5 {print $1}' proteomics_gene_level.tsv \
  > proteomics_genes_ge5_unique_peptides.txt

wc -l proteomics_genes_ge1_unique_peptide.txt
wc -l proteomics_genes_ge2_unique_peptides.txt
wc -l proteomics_genes_ge3_unique_peptides.txt
wc -l proteomics_genes_ge5_unique_peptides.txt
```

Observed proteomics evidence:

```text
>=1 unique peptide: 489 genes
>=2 unique peptides: 205 genes
>=3 unique peptides: 112 genes
>=5 unique peptides: 38 genes
```

---

## 12. Build multi-omics gene evidence matrix

Required input files:

```text
Tt_Final_Annotation.gtf
counts/Tt_DRS_gene_counts.s0.txt
SR_counts/Tt_HVFMLDSX2.s0.txt
LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt
proteomics_gene_level.tsv
```

### 12.1 Check required files

```bash
for f in \
Tt_Final_Annotation.gtf \
counts/Tt_DRS_gene_counts.s0.txt \
SR_counts/Tt_HVFMLDSX2.s0.txt \
LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt \
proteomics_gene_level.tsv
do
  echo "== $f =="
  test -s "$f" && echo "OK" || echo "MISSING"
done
```

### 12.2 Create matrix builder script

```bash
cat > build_multiomics_gene_evidence_matrix.py <<'PY'
#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict


def gene_sort_key(g):
    m = re.search(r"(\d+)", g)
    if m:
        return (re.sub(r"\d+", "", g), int(m.group(1)))
    return (g, 0)


def extract_gene_ids_from_gtf(gtf_path):
    genes = set()
    with open(gtf_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in {"transcript", "gene", "exon", "CDS"}:
                continue
            m = re.search(r'gene_id "([^"]+)"', parts[8])
            if m:
                genes.add(m.group(1))
    return sorted(genes, key=gene_sort_key)


def read_featurecounts_counts(path):
    counts = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#") or row[0] == "Geneid":
                continue
            if len(row) < 7:
                continue
            gene = row[0].strip()
            try:
                count = int(float(row[-1]))
            except ValueError:
                count = 0
            counts[gene] = count
    return counts


def read_proteomics_gene_level(path):
    proteomics = {}
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("GeneID", "").strip()
            if not gid:
                continue
            def to_int(x):
                try:
                    return int(float(str(x).strip()))
                except Exception:
                    return 0
            def to_float(x):
                try:
                    return float(str(x).strip())
                except Exception:
                    return 0.0
            proteomics[gid] = {
                "Max_unique_peptides": to_int(row.get("Max_unique_peptides", 0)),
                "Sum_unique_peptides_reported": to_int(row.get("Sum_unique_peptides_reported", 0)),
                "Proteomics_total_PSMs": to_int(row.get("Total_PSMs", 0)),
                "Proteomics_total_peptides_reported": to_int(row.get("Total_peptides_reported", 0)),
                "Proteomics_max_coverage_percent": to_float(row.get("Max_coverage_percent", 0)),
                "Proteomics_max_Mascot_score": to_float(row.get("Max_Mascot_score", 0)),
                "Best_protein_by_unique_peptides": row.get("Best_protein_by_unique_peptides", "").strip(),
            }
    return proteomics


def yesno(cond):
    return "yes" if cond else "no"


def classify_evidence(drs, sr, lr, prot_unique):
    rna_platforms = sum([drs > 0, sr > 0, lr > 0])
    any_rna = rna_platforms > 0
    prot = prot_unique >= 1
    prot2 = prot_unique >= 2

    if any_rna and prot2:
        if rna_platforms == 3:
            return "three_RNA_and_high_confidence_proteomics_supported"
        return "high_confidence_multiomics_supported"
    if any_rna and prot:
        if rna_platforms == 3:
            return "three_RNA_and_proteomics_supported"
        return "RNA_and_proteomics_supported"
    if prot and not any_rna:
        return "proteomics_only"
    if rna_platforms == 3:
        return "three_RNA_supported"
    if rna_platforms == 2:
        return "multi_RNA_supported"
    if rna_platforms == 1:
        return "single_RNA_supported"
    return "no_detected_omics_support"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--drs", required=True)
    parser.add_argument("--sr", required=True)
    parser.add_argument("--lr", required=True)
    parser.add_argument("--proteomics", required=True)
    parser.add_argument("-o", "--output", default="multiomics_gene_evidence_matrix.tsv")
    args = parser.parse_args()

    genes = extract_gene_ids_from_gtf(args.gtf)
    drs_counts = read_featurecounts_counts(args.drs)
    sr_counts = read_featurecounts_counts(args.sr)
    lr_counts = read_featurecounts_counts(args.lr)
    prot = read_proteomics_gene_level(args.proteomics)

    out_cols = [
        "GeneID", "DRS_count", "DRS_support", "DRS_ge5", "DRS_ge10",
        "SR_cDNA_count", "SR_cDNA_support", "SR_cDNA_ge5", "SR_cDNA_ge10",
        "LR_cDNA_count", "LR_cDNA_support", "LR_cDNA_ge5", "LR_cDNA_ge10",
        "RNA_platforms_supported", "Any_RNA_support", "All_three_RNA_support",
        "Proteomics_support", "Max_unique_peptides", "Sum_unique_peptides_reported",
        "Proteomics_total_PSMs", "Proteomics_total_peptides_reported",
        "Proteomics_max_coverage_percent", "Proteomics_max_Mascot_score",
        "Best_protein_by_unique_peptides",
        "High_confidence_proteomics_2_unique_peptides",
        "Strong_proteomics_3_unique_peptides",
        "Very_strong_proteomics_5_unique_peptides",
        "Multiomics_support", "Evidence_class",
    ]

    summary = defaultdict(int)

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_cols)
        writer.writeheader()

        for gid in genes:
            drs = drs_counts.get(gid, 0)
            sr = sr_counts.get(gid, 0)
            lr = lr_counts.get(gid, 0)
            p = prot.get(gid, {})
            max_unique = int(p.get("Max_unique_peptides", 0))
            prot_support = max_unique >= 1
            prot_ge2 = max_unique >= 2
            rna_platforms = sum([drs > 0, sr > 0, lr > 0])
            any_rna = rna_platforms > 0
            all_three_rna = rna_platforms == 3
            multiomics = any_rna and prot_support
            evidence_class = classify_evidence(drs, sr, lr, max_unique)

            row = {
                "GeneID": gid,
                "DRS_count": drs,
                "DRS_support": yesno(drs > 0),
                "DRS_ge5": yesno(drs >= 5),
                "DRS_ge10": yesno(drs >= 10),
                "SR_cDNA_count": sr,
                "SR_cDNA_support": yesno(sr > 0),
                "SR_cDNA_ge5": yesno(sr >= 5),
                "SR_cDNA_ge10": yesno(sr >= 10),
                "LR_cDNA_count": lr,
                "LR_cDNA_support": yesno(lr > 0),
                "LR_cDNA_ge5": yesno(lr >= 5),
                "LR_cDNA_ge10": yesno(lr >= 10),
                "RNA_platforms_supported": rna_platforms,
                "Any_RNA_support": yesno(any_rna),
                "All_three_RNA_support": yesno(all_three_rna),
                "Proteomics_support": yesno(prot_support),
                "Max_unique_peptides": max_unique,
                "Sum_unique_peptides_reported": p.get("Sum_unique_peptides_reported", 0),
                "Proteomics_total_PSMs": p.get("Proteomics_total_PSMs", 0),
                "Proteomics_total_peptides_reported": p.get("Proteomics_total_peptides_reported", 0),
                "Proteomics_max_coverage_percent": p.get("Proteomics_max_coverage_percent", 0),
                "Proteomics_max_Mascot_score": p.get("Proteomics_max_Mascot_score", 0),
                "Best_protein_by_unique_peptides": p.get("Best_protein_by_unique_peptides", ""),
                "High_confidence_proteomics_2_unique_peptides": yesno(prot_ge2),
                "Strong_proteomics_3_unique_peptides": yesno(max_unique >= 3),
                "Very_strong_proteomics_5_unique_peptides": yesno(max_unique >= 5),
                "Multiomics_support": yesno(multiomics),
                "Evidence_class": evidence_class,
            }
            writer.writerow(row)

            summary["total_genes"] += 1
            summary["DRS_support"] += int(drs > 0)
            summary["SR_cDNA_support"] += int(sr > 0)
            summary["LR_cDNA_support"] += int(lr > 0)
            summary["Any_RNA_support"] += int(any_rna)
            summary["All_three_RNA_support"] += int(all_three_rna)
            summary["Proteomics_support"] += int(prot_support)
            summary["Proteomics_ge2"] += int(prot_ge2)
            summary["RNA_and_proteomics"] += int(any_rna and prot_support)
            summary["Three_RNA_and_proteomics"] += int(all_three_rna and prot_support)
            summary["Three_RNA_and_high_confidence_proteomics"] += int(all_three_rna and prot_ge2)

    print(f"Output written: {args.output}")
    print("Summary:")
    for k in [
        "total_genes", "DRS_support", "SR_cDNA_support", "LR_cDNA_support",
        "Any_RNA_support", "All_three_RNA_support", "Proteomics_support",
        "Proteomics_ge2", "RNA_and_proteomics", "Three_RNA_and_proteomics",
        "Three_RNA_and_high_confidence_proteomics",
    ]:
        print(f"{k}\t{summary[k]}")


if __name__ == "__main__":
    main()
PY

chmod +x build_multiomics_gene_evidence_matrix.py
```

### 12.3 Run matrix builder

```bash
python build_multiomics_gene_evidence_matrix.py \
  --gtf Tt_Final_Annotation.gtf \
  --drs counts/Tt_DRS_gene_counts.s0.txt \
  --sr SR_counts/Tt_HVFMLDSX2.s0.txt \
  --lr LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt \
  --proteomics proteomics_gene_level.tsv \
  -o multiomics_gene_evidence_matrix.tsv

head -5 multiomics_gene_evidence_matrix.tsv
wc -l multiomics_gene_evidence_matrix.tsv
```

Expected matrix size:

```text
33,725 genes + header = 33,726 lines
```

### 12.4 Matrix summary statistics

```bash
awk -F'\t' 'NR>1 {count[$29]++} END{for (c in count) printf "%s\t%d\n", c, count[c]}' \
  multiomics_gene_evidence_matrix.tsv | sort
```

Observed evidence class summary:

```text
multi_RNA_supported                                      1,608
no_detected_omics_support                              11,423
RNA_and_proteomics_supported                                5
single_RNA_supported                                    3,018
three_RNA_and_high_confidence_proteomics_supported        205
three_RNA_and_proteomics_supported                        279
three_RNA_supported                                    17,187
```

Observed global multi-omics summary:

```text
total_genes: 33,725
DRS_support: 18,574
SR_cDNA_support: 19,540
LR_cDNA_support: 21,143
Any_RNA_support: 22,302
All_three_RNA_support: 17,671
Proteomics_support: 489
Proteomics_ge2: 205
RNA_and_proteomics: 489
Three_RNA_and_proteomics: 484
Three_RNA_and_high_confidence_proteomics: 205
```

---

## 13. Figure 6 source data

Figure 6 contains:

- A. DRS / SR cDNA / LR cDNA mapping rate
- B. Gene support by RNA evidence
- C. Multi-omics evidence classes
- D. Proteomics-supported gene evidence classes

Source data:

```text
Figure6_multiomics_support_source_data.tsv
```

Final figure:

```text
Figure6_multiomics_support_Tt_tenax.png
```

---

## 14. Final validation numbers for manuscript

### 14.1 RNA evidence

```text
Total genes: 33,725
DRS-supported genes: 18,574 / 33,725 = 55.1%
SR cDNA-supported genes: 19,540 / 33,725 = 57.9%
LR cDNA-supported genes: 21,143 / 33,725 = 62.7%
Any RNA-supported genes: 22,302 / 33,725 = 66.1%
All three RNA-supported genes: 17,671 / 33,725 = 52.4%
```

### 14.2 Proteomics evidence

```text
Proteomics-supported genes, >=1 unique peptide: 489
High-confidence proteomics-supported genes, >=2 unique peptides: 205
Strong proteomics-supported genes, >=3 unique peptides: 112
Very strong proteomics-supported genes, >=5 unique peptides: 38
RNA + proteomics-supported genes: 489
Three RNA + proteomics-supported genes: 484
Three RNA + high-confidence proteomics-supported genes: 205
```

### 14.3 Suggested Technical Validation paragraph

```text
Multi-omics evidence supported a substantial fraction of the final annotation. Among 33,725 predicted genes, 22,302 (66.1%) were supported by at least one RNA sequencing platform, and 17,671 (52.4%) were supported by all three RNA datasets. Proteomic evidence supported 489 predicted genes, all of which also had RNA support. Of these, 484 genes were supported by all three RNA platforms, including 205 genes with at least two unique peptides. These results provide multi-layer transcript- and protein-level support for the final gene models.
```

---

## 15. Recommended files to archive

```text
# Genome and annotation
ref/draft.polish2.bs6.fasta
Tt_Final_Annotation.gff3
Tt_Final_Annotation.gtf

# Basecalled reads
basecall/Tt_DRS_merged.hac.bam
basecall/Tt_DRS_merged.hac.fastq.gz
basecall/Tt_cDNA_merged.sup.bam
basecall/Tt_cDNA_merged.sup.fastq.gz
basecall/Tt_cDNA_merged.sup.Q10.min200.max50k.fastq.gz

# Alignments
aln/Tt_DRS_merged.hac.genome.bam
SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam
LRcDNA_aln/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.bam

# Counts
counts/Tt_DRS_gene_counts.s0.txt
SR_counts/Tt_HVFMLDSX2.s0.txt
LRcDNA_counts/Tt_cDNA_gene_counts.Q10.min200.max50k.s0.txt

# Proteomics
20260312-PT-15_100-150min-tt.txt
proteomics_gene_level.tsv
proteomics_genes_ge1_unique_peptide.txt
proteomics_genes_ge2_unique_peptides.txt
proteomics_genes_ge3_unique_peptides.txt
proteomics_genes_ge5_unique_peptides.txt

# Multi-omics tables and figure
multiomics_gene_evidence_matrix.tsv
Figure6_multiomics_support_source_data.tsv
Figure6_multiomics_support_Tt_tenax.png
```

---

## 16. Notes and cautions

1. `featureCounts -s 0` was used as the primary mode across DRS, SR cDNA, and LR cDNA because it consistently gave the highest gene assignment or the most conservative annotation-support result.
2. DRS and LR cDNA are long-read transcriptomic evidence, but they should not be interpreted as full transcript isoform resolution without additional transcript reconstruction and curation.
3. The short-read splice junction analysis identified 44 exact annotated intron matches. Because no Sanger sequencing or targeted RT-PCR validation was performed, these are reported as supportive transcriptomic observations only.
4. Proteomics results were interpreted as complementary protein-level support, not as a comprehensive proteome survey.
5. The matrix and Figure 6 should be treated as annotation-support evidence rather than differential expression or quantitative proteomics analysis.
