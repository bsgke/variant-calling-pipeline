#!/bin/bash -ue
trimmomatic PE -phred33 input.1 data/SRR957824_500K_R1.fastq.gz_paired.fq.gz data/SRR957824_500K_R1.fastq.gz_unpaired.fq.gz  data/SRR957824_500K_R1.fastq.gz_paired.fq.gz data/SRR957824_500K_R1.fastq.gz_unpaired.fq.gz ILLUMINACLIP:/home/icipe/miniconda3/envs/variant_calling/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
