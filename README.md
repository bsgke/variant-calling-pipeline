# variant-calling-pipeline
A basic variant calling pipeline 

perform variant calling with the provided script, you can make changes as you see fit.

once it works convert it to nextflow, Good luck

read data:

curl -O -J -L https://osf.io/shqpv/download

curl -O -J -L https://osf.io/9m3ch/download

The dataset you will be working with is from an Illumina MiSeq dataset. The sequenced organism is an
enterohaemorrhagic E. coli (EHEC) of the serotype O157, a potentially fatal gastrointestinal pathogen. The sequenced bacterium was part of an outbreak investigation in the St. Louis area, USA in 2011. The sequencing was done as
paired-end 2x150bp,only a subset of the data is provided with the above links.

reference data: 

curl -O -J -L https://osf.io/rnzbe/download

This reference contains the sequence of the pO157 plasmid from the Sakai outbreak strain of E. coli O157.
