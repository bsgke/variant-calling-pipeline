#making directory
mkdir fastqc 
# fastqc install 
conda install -c bioconda fastqc
# fastqc
fastqc input/*.gz -o input/fastqc/
#multiqc install  
conda install -c bioconda multiqc
# trimmomatics 
conda install -c bioconda trimmomatic
#trimming of reads:Trimmomatic
trimmomatic PE -phred33 input/SRR957824_500K_R1.fastq.gz input/SRR957824_500K_R2.fastq.gz trimmed/SRR957824_500K_R1_paired.fq.gz trimmed/SRR957824_500K_R1_unpaired.fq.gz trimmed/SRR957824_500K_R2_paired.fq.gz trimmed/SRR957824_500K_R2_unpaired.fq.gz ILLUMINACLIP:Truseq3.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#fastqc
fastqc trimmed/*.gz -o trimmed/fastqc/

# bwa installation 
conda install -c bioconda bwa
# indexing reference 
mkdir -p alignment 
bwa  index  reference/pO157_Sakai.fasta.gz 
#alignment: align the trimmed reads with the reference genome
bwa mem reference/pO157_Sakai.fasta.gz trimmed/SRR957824_500K_R1_paired.fq.gz trimmed/SRR957824_500K_R1_paired.fq.gz > alignment/alignment.sam

# sam installation
conda install -c bioconda samtools
#SAM TO BAM
samtools view -Sb alignment/alignment.sam > alignment/alignment.bam
#get the stats for your alignment 
samtools flagstat alignment/alignment.bam

#sambamba installation
conda install -c bioconda sambamba
# removing duplicates 
sambamba markdup -r alignment/alignment.bam alignment/alignment_dedup.bam
#sort bam file
samtools sort alignment/alignment_dedup.bam > alignment/alignment_sorted.bam
#index bam file 
samtools index alignment/alignment_sorted.bam

# VARIANT CALLING STEP 
#making directory
mkdir bcf
#installation bcf 
conda install -c bioconda bcftools
#bcftools for generating the bcf files,ploidy is either haploid or diploid
bcftools mpileup -Ou -A -f reference/pO157_Sakai.fasta.gz alignment/alignment_sorted.bam | bcftools call  -mv -Ob -o bcf/calls.bcf
#bcftools for generating the vcf files 
bcftools view bcf/calls.bcf | vcfutils.pl varFilter -> bcf/bcftools.vcf

#variant calling with freebayes:ploidy is either haploid or diploid
# install freebayes 
conda install -c bioconda freebayes
#variant calling with freebayes:ploidy is either haploid or diploid
freebayes -f reference/pO157_Sakai.fasta.gz --ploidy 2 alignment/alignment_sorted.bam  > freebayes/freebayes.vcf

# check statistics 
#install tabidx 
conda install -c bioconda tabix
#statistics
bgzip freebayes/freebayes.vcf
bgzip bcf/bcftools.vcf 
tabix -p vcf freebayes/freebayes.vcf.gz
tabix -p vcf bcf/bcftools.vcf.gz
bcftools stats -F reference/pO157_Sakai.fasta.gz -s - bcf/bcftools.vcf.gz > bcf/bcf.stats
bcftools stats -F reference/pO157_Sakai.fasta.gz -s - freebayes/freebayes.vcf.gz > freebayes/freebayes.stats


