#indexed the reference files for GATK Analysis
samtools faidx pO157_Sakai.fasta # generates fai indexed files

#GENERATING VCF FILE USING GATK 
#navigate to GATK/home/ckmwangi/Downloads/gatk-4.1.2.0/gatk-4.1.2.0
#Editing the files used for analysis in GATK 
#generates the dict files for indexing
gatk CreateSequenceDictionary --REFERENCE pO157_Sakai.fasta 

 
#Adding Read Groups for the bam files generated
gatk AddOrReplaceReadGroups --INPUT SRR2584857.sorted.bam  --OUTPUT seqcleaned_sorted_groups.bam --RGID 4 --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 

#index the group.bam file:
samtools index seqcleaned_sorted_groups.bam

#call raw variants:
gatk HaplotypeCaller -I seqcleaned_sorted_groups.bam -O raw_gatk_variants.vcf -R pO157_Sakai.fasta

#create SNP & INDEL knownsites if none are available:
gatk SelectVariants -R pO157_Sakai.fasta -V raw_gatk_variants.vcf -select-type SNP -O raw_snps.vcf
gatk SelectVariants -R pO157_Sakai.fasta -V raw_gatk_variants.vcf -select-type INDEL -O raw_indels.vcf

#VariantFiltration for quality variants:
gatk VariantFiltration -R pO157_Sakai.fasta -V raw_snps.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O filtered_snps.vcf

gatk VariantFiltration -R pO157_Sakai.fasta -V raw_indels.vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 2.0' --filter-name "basic_snp_filter" -O filtered_indels.vcf

#base-recalibration:
gatk BaseRecalibrator -R pO157_Sakai.fasta -I seqcleaned_sorted_groups.bam -O Recalibrator_report.grp --known-sites filtered_snps.vcf --known-sites filtered_indels.vcf

gatk ApplyBQSR -bqsr Recalibrator_report.grp -I seqcleaned_sorted_groups.bam -O seqcleaned_sorted_groups_recalibrated.bam

#finally call variants:
gatk HaplotypeCaller -I seqcleaned_sorted_groups_recalibrated.bam -O gatk.vcf -R pO157_Sakai.fasta

#additional analysis:

#bcftools view -v snps bcftools.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -#k1rn

#Annotation:
bgzip gatk.vcf
tabix -p vcf gatk.vcf.gz
bcftools stats -F pO157_Sakai.fasta -s - gatk.vcf.gz > gatk.stats
#use multiqc to see stats including graphs,variant types and numbers
multiqc .

