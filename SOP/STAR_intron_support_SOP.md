# SOP: STAR-supported intron annotation GFF3/GTF

Project: Trichomonas tenax genome annotation validation

Repository note: this repository stores the command SOP, the annotation helper script, and a compact table of supported introns. STAR BAMs, STAR indexes, and other regenerated outputs remain external workflow products.

Purpose: Use Illumina short-read RNA-seq splice junctions from STAR to identify annotated introns with exact coordinate and strand support, then generate STAR-support-enhanced GFF3/GTF files.

Final products:

- `Express/Tt_Final_Annotation.STAR_intron_support.gff3`
- `Express/Tt_Final_Annotation.STAR_intron_support.gtf`
- `Express/Tt_Final_Annotation.STAR_intron_support.summary.tsv`
- `Express/Tt_Final_Annotation.STAR_intron_support.supported_introns.tsv`

Important interpretation:

- The GFF3 keeps all 8,500 annotated intron features and adds STAR support attributes.
- The GTF keeps the original gene/transcript/exon/CDS annotation and writes explicit intron rows only for the 44 STAR-supported introns.
- STAR support means exact match of contig, intron start, intron end, and strand between annotated intron and STAR `SJ.out.tab`.
- The stricter `min5` set means exact STAR-supported introns with total junction reads >= 5.

---

## 0. Working directory and environment

Run all commands from the repository root:

```bash
PROJECT="${PROJECT:-/path/to/Tt_Genome}"
cd "$PROJECT"
```

Recommended tools:

```text
bash
awk
sort
join
wc
STAR
samtools
python3
```

If rebuilding STAR alignment from trimmed reads, activate the RNA-seq environment used for the multi-omics SOP.

```bash
conda activate srna
```

---

## 1. Input files

Required annotation files:

```text
Express/Tt_Final_Annotation.gff3
Express/Tt_Final_Annotation.gtf
```

Required genome FASTA for STAR alignment:

```text
Express/ref/draft.polish2.bs6.fasta
```

Trimmed short-read cDNA FASTQ files:

```text
Express/SR_trim/Tt_HVFMLDSX2.trim.R1.fastq.gz
Express/SR_trim/Tt_HVFMLDSX2.trim.R2.fastq.gz
```

Existing STAR junction file, if STAR alignment was already completed:

```text
Express/SR_aln/Tt_HVFMLDSX2.SJ.out.tab
```

Check inputs:

```bash
ls -lh \
  Express/Tt_Final_Annotation.gff3 \
  Express/Tt_Final_Annotation.gtf \
  Express/ref/draft.polish2.bs6.fasta \
  Express/SR_aln/Tt_HVFMLDSX2.SJ.out.tab
```

---

## 2. Optional: rebuild STAR index

Skip this section if `Express/STAR_index/Tt_genome/` already exists and was generated from the final genome and final GTF.

```bash
mkdir -p Express/STAR_index/Tt_genome

STAR \
  --runThreadN 36 \
  --runMode genomeGenerate \
  --genomeDir Express/STAR_index/Tt_genome \
  --genomeFastaFiles Express/ref/draft.polish2.bs6.fasta \
  --sjdbGTFfile Express/Tt_Final_Annotation.gtf \
  --sjdbOverhang 149
```

---

## 3. Optional: rebuild STAR short-read cDNA alignment

Skip this section if `Express/SR_aln/Tt_HVFMLDSX2.SJ.out.tab` already exists.

STAR BAM sorting can hit temporary-file or open-file limits, so this workflow writes an unsorted BAM first and sorts with `samtools`.

```bash
mkdir -p Express/SR_aln
ulimit -n 65535

rm -rf Express/SR_aln/Tt_HVFMLDSX2._STARtmp

STAR \
  --runThreadN 20 \
  --genomeDir Express/STAR_index/Tt_genome \
  --readFilesIn Express/SR_trim/Tt_HVFMLDSX2.trim.R1.fastq.gz Express/SR_trim/Tt_HVFMLDSX2.trim.R2.fastq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix Express/SR_aln/Tt_HVFMLDSX2. \
  --quantMode GeneCounts

samtools sort -@ 12 -m 2G \
  -o Express/SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam \
  Express/SR_aln/Tt_HVFMLDSX2.Aligned.out.bam

samtools index Express/SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam

samtools flagstat Express/SR_aln/Tt_HVFMLDSX2.Aligned.sortedByCoord.out.bam \
  > Express/SR_aln/Tt_HVFMLDSX2.flagstat.txt

cat Express/SR_aln/Tt_HVFMLDSX2.Log.final.out
cat Express/SR_aln/Tt_HVFMLDSX2.flagstat.txt
```

Observed STAR alignment summary in this project:

```text
Input reads: 37,505,899
Uniquely mapped reads: 33,525,278
Uniquely mapped reads %: 89.39%
Splices total: 212,581
Annotated splices: 198,061
```

---

## 4. Extract annotated introns from GFF3

The final GFF3 contains explicit `intron` features. Extract contig, start, end, strand, transcript ID, and intron ID.

```bash
mkdir -p Express/intron_check

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
}' Express/Tt_Final_Annotation.gff3 | sort -k1,1 -k2,2n > Express/intron_check/annotation_introns.tsv

wc -l Express/intron_check/annotation_introns.tsv
head Express/intron_check/annotation_introns.tsv
```

Expected result:

```text
8,500 annotated introns
```

---

## 5. Convert STAR SJ.out.tab to intron coordinates

STAR `SJ.out.tab` columns used here:

- column 1: contig
- column 2: intron start
- column 3: intron end
- column 4: strand code, where `1=+`, `2=-`, `0=unstranded`
- column 6: annotated junction flag
- column 7: uniquely mapping reads crossing junction
- column 8: multi-mapping reads crossing junction

Convert strand code and compute total junction reads.

```bash
awk '
BEGIN{OFS="\t"}
{
  s="."
  if($4==1) s="+"
  else if($4==2) s="-"
  total=$7+$8
  print $1,$2,$3,s,$7,$8,total,$6
}' Express/SR_aln/Tt_HVFMLDSX2.SJ.out.tab | sort -k1,1 -k2,2n > Express/intron_check/star_junctions.tsv

wc -l Express/intron_check/star_junctions.tsv
head Express/intron_check/star_junctions.tsv
```

Expected result:

```text
4,476 STAR unique junction rows
```

Create the `min5` STAR junction subset.

```bash
awk '$7>=5' Express/intron_check/star_junctions.tsv > Express/intron_check/star_junctions.min5.tsv

wc -l Express/intron_check/star_junctions.min5.tsv
```

Observed result:

```text
2,741 STAR junctions with total reads >= 5
```

---

## 6. Match annotated introns to STAR junctions

Exact support is defined by:

```text
contig:start:end:strand
```

Create keyed tables.

```bash
awk 'BEGIN{OFS="\t"} {key=$1":"$2":"$3":"$4; print key,$0}' \
  Express/intron_check/annotation_introns.tsv > Express/intron_check/annotation_introns.keyed.tsv

awk 'BEGIN{OFS="\t"} {key=$1":"$2":"$3":"$4; print key,$0}' \
  Express/intron_check/star_junctions.tsv > Express/intron_check/star_junctions.keyed.tsv

awk 'BEGIN{OFS="\t"} {key=$1":"$2":"$3":"$4; print key,$0}' \
  Express/intron_check/star_junctions.min5.tsv > Express/intron_check/star_junctions.min5.keyed.tsv
```

Join annotated introns to all exact STAR junctions.

```bash
join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 Express/intron_check/annotation_introns.keyed.tsv) \
  <(sort -k1,1 Express/intron_check/star_junctions.keyed.tsv) \
  > Express/intron_check/annotation_introns_supported_by_STAR.tsv

join -t $'\t' -v 1 -1 1 -2 1 \
  <(sort -k1,1 Express/intron_check/annotation_introns.keyed.tsv) \
  <(sort -k1,1 Express/intron_check/star_junctions.keyed.tsv) \
  > Express/intron_check/annotation_introns_NOT_supported_by_STAR.tsv
```

Join annotated introns to the stricter exact STAR junctions with total reads >= 5.

```bash
join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 Express/intron_check/annotation_introns.keyed.tsv) \
  <(sort -k1,1 Express/intron_check/star_junctions.min5.keyed.tsv) \
  > Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv
```

Check counts.

```bash
wc -l \
  Express/intron_check/annotation_introns.tsv \
  Express/intron_check/star_junctions.tsv \
  Express/intron_check/star_junctions.min5.tsv \
  Express/intron_check/annotation_introns_supported_by_STAR.tsv \
  Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv \
  Express/intron_check/annotation_introns_NOT_supported_by_STAR.tsv
```

Observed result:

```text
8,500 annotated introns
4,476 STAR unique junction rows
2,741 STAR junctions with total reads >= 5
44 annotated introns with exact STAR support
23 annotated introns with exact STAR support and total reads >= 5
8,456 annotated introns without exact STAR support
```

---

## 7. Summarize intron support by transcript

Count total annotated introns per transcript.

```bash
awk 'BEGIN{OFS="\t"} {n[$5]++} END{for(t in n) print t,n[t]}' \
  Express/intron_check/annotation_introns.tsv | sort -k1,1 > Express/intron_check/transcript_total_introns.tsv
```

Count STAR-supported introns per transcript.

```bash
awk 'BEGIN{OFS="\t"} {n[$6]++} END{for(t in n) print t,n[t]}' \
  Express/intron_check/annotation_introns_supported_by_STAR.tsv | sort -k1,1 > Express/intron_check/transcript_supported_introns.tsv
```

Create transcript-level support summary.

```bash
join -t $'\t' -a 1 -e 0 -o '1.1 1.2 2.2' \
  Express/intron_check/transcript_total_introns.tsv \
  Express/intron_check/transcript_supported_introns.tsv \
  > Express/intron_check/transcript_intron_support_summary.tsv
```

Check transcript-level counts.

```bash
wc -l \
  Express/intron_check/transcript_total_introns.tsv \
  Express/intron_check/transcript_supported_introns.tsv \
  Express/intron_check/transcript_intron_support_summary.tsv

awk '$3>0 {n++} END{print "transcripts_with_STAR_supported_introns:", n}' \
  Express/intron_check/transcript_intron_support_summary.tsv

awk '$2==$3 {n++} END{print "transcripts_all_introns_supported:", n}' \
  Express/intron_check/transcript_intron_support_summary.tsv

awk '$3>0 && $2>$3 {n++} END{print "transcripts_partial_intron_support:", n}' \
  Express/intron_check/transcript_intron_support_summary.tsv
```

Observed result:

```text
6,187 intron-containing transcripts
42 transcripts with >=1 STAR-supported intron
34 transcripts with all annotated introns supported
8 transcripts with partial intron support
```

---

## 8. Add STAR support attributes to GFF3/GTF

The script used for the final annotation-support files is:

```text
scripts/add_star_intron_support_to_annotations.py
```

Run with default paths.

```bash
python3 scripts/add_star_intron_support_to_annotations.py
```

Equivalent explicit command:

```bash
python3 scripts/add_star_intron_support_to_annotations.py \
  --gff3 Express/Tt_Final_Annotation.gff3 \
  --gtf Express/Tt_Final_Annotation.gtf \
  --supported Express/intron_check/annotation_introns_supported_by_STAR.tsv \
  --supported-min5 Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv \
  --out-prefix Express/Tt_Final_Annotation.STAR_intron_support
```

Output files:

```text
Express/Tt_Final_Annotation.STAR_intron_support.gff3
Express/Tt_Final_Annotation.STAR_intron_support.gtf
Express/Tt_Final_Annotation.STAR_intron_support.summary.tsv
Express/Tt_Final_Annotation.STAR_intron_support.supported_introns.tsv
```

GFF3 attributes added to gene and mRNA features when the model contains annotated introns:

```text
STAR_annotated_introns
STAR_supported_introns
STAR_supported_introns_min5
STAR_any_intron_support
STAR_all_introns_supported
```

GFF3 attributes added to intron features:

```text
STAR_junction_support
STAR_unique_reads
STAR_multimap_reads
STAR_total_reads
STAR_annotated_junction
STAR_min5_support
```

GTF attributes added to transcript features:

```text
star_annotated_introns
star_supported_introns
star_supported_introns_min5
star_any_intron_support
star_all_introns_supported
```

GTF explicit intron rows are written only for STAR-supported introns and include:

```text
gene_id
transcript_id
intron_id
star_junction_support
star_unique_reads
star_multimap_reads
star_total_reads
star_annotated_junction
star_min5_support
```

---

## 9. Validate final output

Check script syntax.

```bash
python3 -m py_compile scripts/add_star_intron_support_to_annotations.py
```

Check output files.

```bash
ls -lh \
  Express/Tt_Final_Annotation.STAR_intron_support.gff3 \
  Express/Tt_Final_Annotation.STAR_intron_support.gtf \
  Express/Tt_Final_Annotation.STAR_intron_support.summary.tsv \
  Express/Tt_Final_Annotation.STAR_intron_support.supported_introns.tsv
```

Check STAR support counts in GFF3.

```bash
rg -c 'STAR_junction_support=exact' Express/Tt_Final_Annotation.STAR_intron_support.gff3
rg -c 'STAR_min5_support=yes' Express/Tt_Final_Annotation.STAR_intron_support.gff3
```

Expected result:

```text
44
23
```

Check STAR-supported explicit intron rows in GTF.

```bash
rg -c $'\tintron\t' Express/Tt_Final_Annotation.STAR_intron_support.gtf
rg -c 'star_junction_support "exact"' Express/Tt_Final_Annotation.STAR_intron_support.gtf
rg -c 'star_junction_support "no"' Express/Tt_Final_Annotation.STAR_intron_support.gtf
```

Expected result:

```text
44
44
0
```

Check final summary.

```bash
cat Express/Tt_Final_Annotation.STAR_intron_support.summary.tsv
```

Expected summary:

```text
annotated_introns_in_gff3: 8500
star_exact_supported_introns: 44
star_exact_supported_introns_min5: 23
explicit_introns_written_to_gtf: 44
intron_containing_transcripts: 6187
transcripts_with_star_supported_introns: 42
transcripts_all_introns_supported: 34
intron_containing_genes: 6149
genes_with_star_supported_introns: 42
genes_all_introns_supported: 34
```

Inspect example records.

```bash
rg -n -m 5 'STAR_junction_support=exact' Express/Tt_Final_Annotation.STAR_intron_support.gff3
rg -n -m 5 'star_junction_support "exact"' Express/Tt_Final_Annotation.STAR_intron_support.gtf
```

Example supported intron in the final GTF:

```text
g496.t1.intron1
contig_100:1162615-1162687:+
STAR unique reads: 6
STAR total reads: 6
STAR min5 support: yes
```

---

## 10. Manuscript/report wording

Suggested concise wording:

```text
Short-read RNA-seq splice junctions were compared against annotated introns by exact coordinate and strand matching. Among 8,500 annotated introns, 44 introns from 42 transcript models were supported by exact STAR junctions, and 23 had total junction read support >=5. These STAR-supported introns were retained as transcriptomic support evidence. The STAR-supported GTF contains explicit intron rows only for the 44 supported introns, whereas the GFF3 retains all annotated introns with STAR support flags.
```

---

## 11. Output retention

Keep these files with the analysis archive:

```text
Express/intron_check/annotation_introns.tsv
Express/intron_check/star_junctions.tsv
Express/intron_check/star_junctions.min5.tsv
Express/intron_check/annotation_introns_supported_by_STAR.tsv
Express/intron_check/annotation_introns_supported_by_STAR.min5.tsv
Express/intron_check/annotation_introns_NOT_supported_by_STAR.tsv
Express/intron_check/transcript_intron_support_summary.tsv
scripts/add_star_intron_support_to_annotations.py
Express/Tt_Final_Annotation.STAR_intron_support.gff3
Express/Tt_Final_Annotation.STAR_intron_support.gtf
Express/Tt_Final_Annotation.STAR_intron_support.summary.tsv
Express/Tt_Final_Annotation.STAR_intron_support.supported_introns.tsv
```
