#!/bin/bash -ue
trimmomatic PE -phred33 SRR957824_500K_R1.fastq.gz SRR957824_500K_R2.fastq.gz SRR957824_500K_R1_paired.fq.gz SRR957824_500K_R1_unpaired.fq.gz  SRR957824_500K_R1_paired.fq.gz SRR957824_500K_R2_unpaired.fq.gz ILLUMINACLIP:/home/icipe/miniconda3/envs/variant_calling/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
echo SRR957824_500K_R1_fastqc.html SRR957824_500K_R1_fastqc.zip SRR957824_500K_R2_fastqc.html SRR957824_500K_R2_fastqc.zip
