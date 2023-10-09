#!/bin/bash

outdir="results"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
bamlist="muskox_bamlist.txt"
scaffolds="contigs1MB_remap_autosomes_headers.txt"

## Generate GL in beagle format
angsd -b $bamlist -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -P 30 -minInd 54 -minMaf 0.05 \
-rf $scaffolds -out $outdir/muskox_remap2

python pcangsd.py -beagle muskox_remap2.beagle.gz -kinship -sites_save -inbreedSites -o $outdir/muskox_whitelist -threads 20 > $outdir/muskox.log
