#!/bin/bash

outdir="results"
bamlist="muskox_bamlist.txt"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
scaffolds="whitelist.scaffolds"
regions="whitelist.sites"

## Generate GL in beagle format
angsd -b $bamlist -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -P 4 -minInd 54 -minMaf 0.05 \
-rf $scaffolds -sites $regions -out $outdir/muskox_muskoxmap_whitelist
