# RNA_Seq 
#### RNA-seq analysis pipeline
#### 1. Download reference genome and index fasta file hg38p12
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
#### 2. Align reads to reference genome
```
bwa mem -t 4 hg38.p12.fa 253_POS_1_S1_R1_001.fastq.gz 253_POS_1_S1_R2_001.fastq.gz | samtools view -b > POS_1.raw.bam
```
#### 3. Sort and index bam file
```
sambamba sort -t 1 -m 32G -p --tmpdir /scratch -o POS_1.bam POS_1.raw.bam
```
#### 4. Download annotation gtf hg38 gtf
```
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
```
#### 5. Count reads
```
featureCounts.img featureCounts -p -B -C -t exon -g gene_id -a /data/annotation/Homo_sapiens.GRCh38.95.gtf -F GTF -o /data/countsLonardo.txt /data/raw/bam/CTR.1.bam /data/raw/bam/CTR.2.bam /data/raw/bam/CTR.3.bam /data/raw/bam/NODAL.1.bam /data/raw/bam/NODAL.2.bam /data/raw/bam/NODAL.3.bam /data/raw/bam/POS.1.bam /data/raw/bam/POS.2.bam /data/raw/bam/POS.3.bam /data/raw/bam/NEG.1.bam /data/raw/bam/NEG.2.bam /data/raw/bam/NEG.3.bam
```
#### 6. Reorganize file counts.txt (remove first line, change sample names, take only information about reads count)
```
cat countsLonardo.txt | grep -v "^#" >allcounts.txt
```
```
cat allcounts.txt  | awk '{print $1, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18} ' | tr " " "\t"  > Lonardo_counts.txt
```
#### 7. Create samples and contrast files
```
samples.tsv
CTR	CTR.1
CTR	CTR.2
CTR	CTR.3
NODAL	NODAL.1
NODAL	NODAL.2
NODAL	NODAL.3
POS	POS.1
POS	POS.2
POS	POS.3
NEG	NEG.1
NEG	NEG.2
NEG	NEG.3
```
```
contrast.tsv
CTR	NODAL
POS	NEG
```
#### 8. Run edgeR
```
perl run_DE_analysis.pl  --matrix Lonardo_counts.txt  --method edgeR --dispersion 0.1 --samples_file samples.tsv --contrasts contrast.tsv  --output results
```
#### 9. Run gsea
#### 10. Post edgeR
