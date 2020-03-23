#quality check with fastqc

fastqc *.gz

#trimming of reads:Trimmomatic

trimmomatic PE -phred33 read-R1.fastq.gz read-R2.fastq.gz read_forward_paired.fq.gz read_forward_unpaired.fq.gz read_reverse_paired.fq.gz read_reverse_unpaired.fq.gz ILLUMINACLIP:Truseq3.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#OR:SolexaQA

SolexaQA++ dynamictrim read-R1.fastq.gz read-R2.fastq.gz 

#OR:Sickle

sickle pe  -f read-R1.fastq.gz -r read-R2.fastq.gz -t sanger -o Dir/trimmed_PE_forward.fastq -p sickle/trimmed_PE_reverse.fastq -s trimmed_singles.fastq

#index the reference:
bwa index -p Ref Ref.fasta

#alignment: align the trimmed reads with the reference genome
bwa mem Ref read-R1.trim.fastq read-R2.trim.fastq > bwa.sam

#SAM TO BAM
samtools view -Sb bwa.sam > bwa.bam

#get the stats for your alignment 
samtools flagstat bwa.bam 

#mark and remove duplicates
sambamba markdup -r bwa.bam dedup_bwa.bam


#sort bam file
samtools sort dedup_bwa.bam > sorted_dedup_bwa.bam


#index bam file 
samtools index sorted_dedup_bwa.bam 

#bcftools for generating the bcf files,ploidy is either haploid or diploid
bcftools mpileup -f Ref.fasta sorted_dedup_bwa.bam | bcftools call --ploidy 2 -mv -Ob -o calls.bcf

#bcftools for generating the vcf files 
bcftools view calls.bcf | vcfutils.pl varFilter -> bcftools.vcf

#variant calling with freebayes:ploidy is either haploid or diploid
freebayes -f Ref.fasta --ploidy 2 seqcleaned_sorted.bam  > freebayes.vcf

bgzip bcftools.vcf 

bgzip freebayes.vcf 

tabix -p vcf bcftools.vcf.gz

tabix -p vcf freebayes.vcf.gz

bcftools stats -F Ref.fasta -s - bcftools.vcf.gz > bcf.stats

bcftools stats -F Ref.fasta -s - freebayes.vcf.gz > freebayes.stats

#use multiqc to generate a report for entire pipeline:
multiqc .







