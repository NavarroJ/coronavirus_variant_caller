#!/bin/bash

bwa index NC_045512.2.fasta
bwa mem NC_045512.2.fasta ncov_global.fasta > ncov_global.sam
samtools fixmate -O bam ncov_global.sam ncov_global.bam
samtools sort -O bam -o ncov_global.sorted.bam ncov_global.bam
samtools index ncov_global.sorted.bam
samtools mpileup -f NC_045512.2.fasta -B -R -aa ncov_global.sorted.bam -o ncov_global.mpileup
