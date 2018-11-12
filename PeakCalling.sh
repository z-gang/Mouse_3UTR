#!/bin/bash

#BSUB -J ReadsMap
#BSUB -n 20
#BSUB -R span[ptile=10]
#BSUB -R rusage[mem=8]
#BSUB -W 12:00
#BSUB -o Mapping.%J.stdout
#BSUB -eo Mapping.%J.stderr

module load samtools/1.7

INDEXDIR="/home/zheng/RefGenome/mm10/GenomeHisat2Index"
FA="mm10.fa"
RAWDIR="/data/mayrc/data/zheng/Mathieu/FLAMAND_5168_181005B1"
FILTERDIR="/data/mayrc/data/zheng/Mathieu/Filtered"
CLEANDIR="/data/mayrc/data/zheng/Mathieu/Clean"
BAMDIR="/data/mayrc/zheng/Mathieu/Bam"

cd $LS_SUBCWD

for i in C57__mNH-12_13__S6 C57__mNH-9__S5 C57__Cre__1__mNH-10__S3 C57__Cre__2__mNH-14_15__S4 Mettl3_F_F___Cre__1__mNH-7__S1 Mettl3_F_F___Cre__1__mNH-8__S2
do
###trim reads by quality
trimmomatic SE -threads 20 $RAWDIR/${i}_L001_R1_001.fastq $FILTERDIR/${i}.filtered.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:21
### trim adapters from reads and trim polyA out
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 21 $FILTERDIR/${i}.filtered.fastq.gz | cutadapt -g "A{100}" -m 21 - | cutadapt -g AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --discard-trimmed --overlap=15 - -o $CLEANDIR/${i}.fastq.gz
### map high quality reads to the genome using hisat2
hisat2 -p 20 -x $INDEXDIR/$FA -U $CLEANDIR/${i}.fastq.gz | samtools view -@ 20 -bS - | samtools sort -@ 20 - > $BAMDIR/${i}.sorted.bam
### index bam file
samtools index $BAMDIR/${i}.sorted.bam -@ 20
done

cd /data/mayrc/zheng/Mathieu/Bam

### merge Bam files for each sample
samtools merge C57.sorted.bam C57__mNH-12_13__S6.sorted.bam C57__mNH-9__S5.sorted.bam -@ 6
samtools merge C57_Cre.sorted.bam C57__Cre__1__mNH-10__S3.sorted.bam C57__Cre__2__mNH-14_15__S4.sorted.bam -@ 6
samtools merge Mettle3_Cre.sorted.bam Mettl3_F_F___Cre__1__mNH-7__S1.sorted.bam Mettl3_F_F___Cre__1__mNH-8__S2.sorted.bam -@ 6

### filter bam file to keep only reads with quality score higher than 15
samtools view -h -bq 15 -@ 6 C57.sorted.bam > C57.filtered.bam
samtools view -h -bq 15 -@ 6 C57_Cre.sorted.bam > C57_Cre.filtered.bam
samtools view -h -bq 15 -@ 6 Mettle3_Cre.sorted.bam > Mettle3_Cre.filtered.bam

### index filtered bam files
samtools index C57.filtered.bam -@ 6
samtools index C57_Cre.filtered.bam -@ 6
samtools index Mettle3_Cre.filtered.bam -@ 6

cd $LS_SUBCWD

### Calling peaks
Rscript PeakCalling.mouse.r
