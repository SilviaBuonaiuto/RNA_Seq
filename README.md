# RNA_Seq 
#### RNA-seq analysis pipeline
#### Required tools and programming language
bwa (http://bio-bwa.sourceforge.net/bwa.shtml)
samtools (http://www.htslib.org/)
sambamba (https://lomereiter.github.io/sambamba/)
featureCounts (http://subread.sourceforge.net/)
edgeR(https://bioconductor.org/packages/release/bioc/html/edgeR.html)
R
#### 1. Download reference genome and index fasta file hg38p12
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
#### 2. Align paired ends reads to reference genome and write raw bam file
```
for i in ctrl1 ctrl2 ctrl3 sample1 sample2 sample3; do bwa mem -t 4 hg38.p12.fa $i_R1_001.fastq.gz $i_R2_001.fastq.gz | samtools view -b > $i.raw.bam ; done
```
#### 3. Sort and index bam file
```
for i in ctrl1 ctrl2 ctrl3 sample1 sample2 sample3; do sambamba sort -t 1 -m 32G -p --tmpdir /scratch -o $i.bam $i.raw.bam
```
#### 4. Download annotation gtf hg38 gtf
```
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
```
#### 5. Count reads
```
featureCounts -p -B -C -t exon -g gene_id -a /data/annotation/Homo_sapiens.GRCh38.95.gtf -F GTF -o readCounts.txt ctrl1.bam ctrl2.bam ctrl3.bam sample1.bam sample2.bam sample3.bam
```
#### 6. Reorganize file counts.txt (remove first line, change sample names, take only information about reads count)
```
cat readCounts.txt | grep -v "^#" > allcounts.txt
```
```
cat allcounts.txt  | awk '{print $1, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18} ' | tr " " "\t"  > finalReadCounts.txt
```
#### 7. Create samples and contrast files
```
samples.tsv
CTR	ctrl1
CTR	ctrl2
CTR	ctrl3
SAMPLE	sample1
SAMPLE	sample2
SAMPLE	sample3
```
```
contrast.tsv
CTR	SAMPLE
```
#### 8. Run edgeR
```
perl run_DE_analysis.pl  --matrix finalReadCounts.txt --method edgeR --dispersion 0.1 --samples_file samples.tsv --contrasts contrast.tsv  --output results
```
#### 9. Run gsea
#### 10. Post edgeR