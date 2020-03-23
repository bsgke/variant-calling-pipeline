#!/usr/bin/env nextflow

// this is a script design in our study group
/*
 * Define the default parameters 
 */
params.inputs = "data/*_R{1,2}.fastq.gz"
params.ref = "reference/*.fasta.gz"
Channel
      .fromFilePairs(params.inputs)
      .into{input1_ch;input2_ch} 

REF= file(params.ref)

process fastqc{
        publishDir 'data/fastqc'
        input:
        set replicate_id, file(reads) from input1_ch
        
        output:
        file ("${replicate_id}*.*") into fastqc_let
        
        script:
        "fastqc $reads" 
}

process trim{
        publishDir 'trimmed'
        input:
        set pair_id, file(input1) from input2_ch
        file clean from fastqc_let
        output:
	file ("${pair_id}*paired.fq.gz") into paired_files
	file ("${pair_id}*unpaired.fq.gz") into unpaired_files
    
        script:
        """
         trimmomatic PE -phred33 $input1 ${pair_id}_R1_paired.fq.gz ${pair_id}_R1_unpaired.fq.gz  ${pair_id}_R2_paired.fq.gz ${pair_id}_R2_unpaired.fq.gz ILLUMINACLIP:/home/icipe/miniconda3/envs/variant_calling/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
         echo $clean
       """

}

process clean_fastqc{
        tag "$pair_id"
        publishDir 'trimmed/fastqc'
        input:
        set pair_id, file(input3) from paired_files
        output:
        file ("*") into clean_ch
        file "input3" into orig_ch1
        script:
        
       "fastqc $input3"
}

process index{
        publishDir 'reference'

        input:
        file ref from REF
        file clean2 from clean_ch
        
        output:
        file ("${ref.baseName}.gz.{amb,ann,bwt,pac,sa}")  into seq_ch
        
        script:
       """
        bwa index $ref
        echo $clean2
        """
}


process alignment{
	tag "$pair_id"
        publishDir 'input/alignment'
        input:
        file reference from REF
        file read_paired from orig_ch1
        output:
        set replicate_id, file ('sam') into output
        script:

       "bwa mem $reference $read_paired > alignment.sam"
}


output.subscribe { print "Done!" }




