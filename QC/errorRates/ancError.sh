#!/bin/bash

output="results"
HD3bam="MQ30bam/HD-3.MQ30.bam"
ref="HD-3.fa.gz"
bamlist="muskox_sheepmap_bamlist.txt"
anc="GCF_002742125.1_Oar_rambouillet_v1.0_genomic_simpleheader.fna"
scaffolds="GCF_002742125.1_Oar_rambouillet_v1.0_genomic_whitelist_autosomes.headers"

angsd -i $HD3bam -doCounts 1 -doFasta 1 -explode 1 -minMapQ 30 -minQ 30 -out $output/HD-3

samtools faidx $ref

angsd -doAncError 1 -anc $anc -ref $ref -nThreads 20 -out $output/muskox_sheepmap_error \
-bam $bamlist -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -rf $scaffolds 
