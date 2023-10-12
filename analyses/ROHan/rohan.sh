#!/bin/bash

output="results"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
contigs="whitelist.scaffolds"
sites="whitelist.sites"
rohan="rohan"
bamfolder="mapping2/results/bam"
outbam="mapping2/results/bam/filt"
samtools="/samtools"

#do not filter for mapping quality
for file in $bamfolder/*.bam
do
	$samtools view -hb -F 3844 -@ 20 $file -o $outbam/$(basename "${file%.*}")_filt.bam
	$samtools index $outbam/$(basename "${file%.*}")_filt.bam
	$rohan --rohmu 2e-5 --map $sites -t 40 -o $(basename "${file%.*}")_rohan $ref $outbam/$(basename "${file%.*}")_filt.bam
done
