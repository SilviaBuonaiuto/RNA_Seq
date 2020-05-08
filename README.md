# RNA_Seq 
#### RNA-seq analysis pipeline
#### 1. Download reference genome and index fasta file hg38p12
#### 2. Align reads to reference genome
#### 3. Sort and index bam file
#### 4. Download annotation gtf hg38 gtf
#### 5. Count reads
#### 6. Reorganize file counts.txt (remove first line, change sample names, take only information about reads count)
#### 7. Create samples and contrast files
#### 8. Run edgeR
#### 9. Run gsea
#### 10. Post edgeR
