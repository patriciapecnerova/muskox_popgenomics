#!/bin/bash

outdir="results"
bamlist="muskox_groups.txt" ## one-column list of links to per-population files that contain a list of bam files in each pop
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
scaffolds="whitelist.scaffolds"
regions="whitelist.sites"

cat $bamlist | while read line
do
angsd -bam "$line" -out $outdir/$(basename "${line%.*}") \
-minMapQ 30 -minQ 30 -uniqueOnly 1 -remove_bads 1 -dosaf 1 -sites $regions -rf $scaffolds -anc $ref -GL 2 -P 2
done

angsd/misc/realSFS -P 10 muskox_bamlist_CaIN.saf.idx > muskox_bamlist_CaIN.sfs 2> muskox_bamlist_CaIN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaIS.saf.idx > muskox_bamlist_CaIS.sfs 2> muskox_bamlist_CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaMW.saf.idx > muskox_bamlist_CaMW.sfs 2> muskox_bamlist_CaMW.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaME.saf.idx > muskox_bamlist_CaME.sfs 2> muskox_bamlist_CaME.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrEN.saf.idx > muskox_bamlist_GrEN.sfs 2> muskox_bamlist_GrEN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx > muskox_bamlist_GrES.sfs 2> muskox_bamlist_GrES.sfs.log

angsd/misc/realSFS -P 10 muskox_bamlist_CaIN.saf.idx muskox_bamlist_CaIS.saf.idx > muskox_bamlist_CaIN.CaIS.sfs 2> muskox_bamlist_CaIN.CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaMW.saf.idx muskox_bamlist_CaIS.saf.idx > muskox_bamlist_CaMW.CaIS.sfs 2> muskox_bamlist_CaMW.CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaME.saf.idx muskox_bamlist_CaIS.saf.idx > muskox_bamlist_CaME.CaIS.sfs 2> muskox_bamlist_CaME.CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrEN.saf.idx muskox_bamlist_CaIS.saf.idx > muskox_bamlist_GrEN.CaIS.sfs 2> muskox_bamlist_GrEN.CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx muskox_bamlist_CaIS.saf.idx > muskox_bamlist_GrES.CaIS.sfs 2> muskox_bamlist_GrES.CaIS.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaMW.saf.idx muskox_bamlist_CaIN.saf.idx > muskox_bamlist_CaMW.CaIN.sfs 2> muskox_bamlist_CaMW.CaIN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaME.saf.idx muskox_bamlist_CaIN.saf.idx > muskox_bamlist_CaME.CaIN.sfs 2> muskox_bamlist_CaME.CaIN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrEN.saf.idx muskox_bamlist_CaIN.saf.idx > muskox_bamlist_GrEN.CaIN.sfs 2> muskox_bamlist_GrEN.CaIN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx muskox_bamlist_CaIN.saf.idx > muskox_bamlist_GrES.CaIN.sfs 2> muskox_bamlist_GrES.CaIN.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_CaME.saf.idx muskox_bamlist_CaMW.saf.idx > muskox_bamlist_CaME.CaMW.sfs 2> muskox_bamlist_CaME.CaMW.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrEN.saf.idx muskox_bamlist_CaMW.saf.idx > muskox_bamlist_GrEN.CaMW.sfs 2> muskox_bamlist_GrEN.CaMW.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx muskox_bamlist_CaMW.saf.idx > muskox_bamlist_GrES.CaMW.sfs 2> muskox_bamlist_GrES.CaMW.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrEN.saf.idx muskox_bamlist_CaME.saf.idx > muskox_bamlist_GrEN.CaME.sfs 2> muskox_bamlist_GrEN.CaME.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx muskox_bamlist_CaME.saf.idx > muskox_bamlist_GrES.CaME.sfs 2> muskox_bamlist_GrES.CaME.sfs.log
angsd/misc/realSFS -P 10 muskox_bamlist_GrES.saf.idx muskox_bamlist_GrEN.saf.idx > muskox_bamlist_GrES.GrEN.sfs 2> muskox_bamlist_GrES.GrEN.sfs.log
