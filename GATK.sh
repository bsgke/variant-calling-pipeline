#!/bin/sh

bwa index pO157_Sakai.fasta

for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
bwa mem pO157_Sakai.fasta  ${infile} ${base}_2.fastq.gz > ${base}.sam
done
	#for infile in  *_R1.fastq 
	#do 
	#base=$(basename ${infile} _R1.fastq) 
	#bwa mem contigs.fa ${infile} ${base}_R2.fastq > ${base}.sam
	#done



for infile in *.sam
do
base=$(basename ${infile} .sam)
samtools view -Sb ${infile} > ${base}.bam
done

parallel samtools flagstat ::: *.bam > flagstats

for infile in *.bam
do
base=$(basename ${infile} .bam)
samtools view -h -b -q 20 ${infile}  > ${base}.q20.bam 
done

for infile in *.q20.bam
do
base=$(basename ${infile} .q20.bam)
sambamba markdup ${infile}  ${base}dedup_q20.bam
done

for infile in *_q20.bam
do
base=$(basename ${infile} _q20.bam)
lofreq viterbi -f pO157_Sakai.fasta ${infile} -o ${base}_realign.bam
done

for infile in *_realign.bam
do
base=$(basename ${infile} _realign.bam)
samtools sort ${infile} ${base}_sorted
done 


samtools faidx pO157_Sakai.fasta

gatk CreateSequenceDictionary --REFERENCE pO157_Sakai.fasta

for infile in *_sorted.bam
do
base=$(basename ${infile} _sorted.bam)
gatk AddOrReplaceReadGroups --INPUT ${infile}  --OUTPUT ${base}_groups.bam --RGID 4 --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 
done

parallel samtools index ::: *_groups.bam


for infile in *_groups.bam
do
base=$(basename ${infile} _groups.bam)
gatk HaplotypeCaller -R pO157_Sakai.fasta -I ${infile} -O ${base}_rawgatk.vcf
done


for infile in *_rawgatk.vcf
do
base=$(basename ${infile} _rawgatk.vcf)
gatk SelectVariants -R pO157_Sakai.fasta -V ${infile} -select-type SNP -O ${base}raw_snps.vcf
gatk SelectVariants -R pO157_Sakai.fasta -V ${infile} -select-type INDEL -O ${base}raw_indels.vcf
done


for infile in *_snps.vcf
do
base=$(basename ${infile} _snps.vcf)
gatk VariantFiltration -R pO157_Sakai.fasta -V ${infile} --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O ${base}filtered_snps.vcf
done

for infile in *_indels.vcf
do
base=$(basename ${infile} _indels.vcf)
gatk VariantFiltration -R pO157_Sakai.fasta -V ${infile} --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 2.0' --filter-name "basic_snp_filter" -O ${base}filtered_indels.vcf
done

for infile in *_groups.bam
do
base=$(basename ${infile} _groups.bam)
gatk BaseRecalibrator -R pO157_Sakai.fasta -I ${infile} -O ${base}Recalibrator_report.grp --known-sites ${base}rawfiltered_snps.vcf --known-sites ${base}rawfiltered_indels.vcf
done

for infile in *Recalibrator_report.grp
do
base=$(basename ${infile} Recalibrator_report.grp)
gatk ApplyBQSR -bqsr ${infile} -I ${base}_groups.bam -O ${base}_recalibrated.bam
done

for infile in *_recalibrated.bam
do
base=$(basename ${infile} _recalibrated.bam)
gatk HaplotypeCaller -I ${infile} -O ${base}_gatk.vcf -R pO157_Sakai.fasta
done

parallel bgzip ::: *_gatk.vcf
parallel tabix -p vcf ::: *_gatk.vcf.gz

for infile in *_gatk.vcf.gz
do
base=$(basename ${infile} _gatk.vcf.gz)
bcftools stats -F pO157_Sakai.fasta -s - ${infile} > ${base}_filtered.stats
done

parallel rtg vcfstats ::: *_gatk.vcf.gz > gatk_vcf.txt

multiqc .


mkdir plots

plot-vcfstats -p plots/ *_filtered.stats
