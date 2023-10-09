#!/bin/bash

output="results"
bamlist="muskox_bamlist.txt"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
contigs="whitelist.scaffolds"
sites="whitelist.sites"

cat ${bamlist} | while read line
do
angsd -i "$line" -anc $ref -dosaf 1 -minMapQ 30 -minQ 30 -P 10 -out $output/$(basename "${line%.*}") \
-uniqueOnly 1 -remove_bads 1 -rf $contigs -sites $sites \
-gl 2 -doMajorMinor 1 -doCounts 1 -setMinDepth 3
angsd/misc/realSFS -fold 1 $output/$(basename "${line%.*}").saf.idx > $output/$(basename "${line%.*}")_est.ml
done

mkdir all

awk '{print FILENAME, $0}' *_est.ml > all/all_est.ml
sed -i 's/_est.ml//g' all/all_est.ml
