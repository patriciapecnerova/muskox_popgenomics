#!/bin/bash

OUTDIR="results"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
whitelistScaffolds="=whitelist.scaffolds"
whitelistSites="whitelist.sites"

######## CaMW
## Generate GL in beagle format
angsd -b muskox_bamlist_CaMW.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_CaMW

zcat $OUTDIR/muskox_CaMW.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 15 -z muskox_bamlist_CaMW.txt -G $OUTDIR/muskox_CaMW.beagle.gz -f freq -O $OUTDIR/muskox_CaMW_ngsrelate

######## CaME
angsd -b muskox_bamlist_CaME.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_CaME

zcat $OUTDIR/muskox_CaME.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 22 -z muskox_bamlist_CaME.txt -G $OUTDIR/muskox_CaME.beagle.gz -f freq -O $OUTDIR/muskox_CaME_ngsrelate

######## CaIS
angsd -b muskox_bamlist_CaIS.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_CaIS

zcat $OUTDIR/muskox_CaIS.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 16 -z muskox_bamlist_CaIS.txt -G $OUTDIR/muskox_CaIS.beagle.gz -f freq -O $OUTDIR/muskox_CaIS_ngsrelate

######## CaIN
angsd -b muskox_bamlist_CaIN.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_CaIN

zcat $OUTDIR/muskox_CaIN.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 27 -z muskox_bamlist_CaIN.txt -G $OUTDIR/muskox_CaIN.beagle.gz -f freq -O $OUTDIR/muskox_CaIN_ngsrelate

######## GrEN
angsd -b muskox_bamlist_GrEN.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_GrEN

zcat $OUTDIR/muskox_GrEN.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 20 -z muskox_bamlist_GrEN.txt -G $OUTDIR/muskox_GrEN.beagle.gz -f freq -O $OUTDIR/muskox_GrEN_ngsrelate


######## GrES
angsd -b muskox_bamlist_GrES.txt -ref $ref -uniqueOnly 1 -remove_bads 1 \
-gl 2 -doGlf 2 -doMajorMinor 1 -minMapQ 30 -minQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 \
-rf $whitelistScaffolds -sites $whitelistSites -out $OUTDIR/muskox_GrES

zcat $OUTDIR/muskox_GrES.mafs.gz | cut -f6 |sed 1d > freq

ngsRelate -p 8 -n 7 -z muskox_bamlist_GrES.txt -G $OUTDIR/muskox_GrES.beagle.gz -f freq -O $OUTDIR/muskox_GrES_ngsrelate
