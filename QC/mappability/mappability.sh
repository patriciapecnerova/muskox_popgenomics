#!/bin/bash

## add links to the software and locations to the files
output="results"
ref="muskox_pseudohap.gapcloser.v1.lenchecked_simpleheaders.fasta"
bamlist="muskox_bamlist.txt"
GEM="GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin"
bedmap="bedops/applications/bed/bedmap/bin/bedmap-typical"
bedtools="bedtools2/bin/bedtools"

##choose your species
species=muskox

## choose kmer size based on the read length of your sequencing reads
kmer=150

## following the pipeline described here: https://evodify.com/gem-mappability/
## start by adding programs to path, because GEM is super finicky
export PATH=$PATH:GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
export PATH=$PATH:GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-indexer_fasta2meta+cont

## indexing the reference
${GEM}/gem-indexer -T 30 -i ${ref} -o ${output}/"${ref##*/}"_gem_index

## generating mappability scores
## allowing for 2 mismatches (1.34% of 150bp read length) and 2 differences
${GEM}/gem-mappability -T 30 -m 0.0134 -e 0.0134 -I ${output}/"${ref##*/}"_gem_index.gem -l 150 -o ${species}_mappability_${kmer}

## converting the GEM file into a bed file
## using tools downloaded here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
## and the last step here: https://github.com/xuefzhao/Reference.Mappability/blob/master/Scripts/bedGraphTobed
## remember to make them executable
${GEM}/gem-2-wig -I ${output}/"${ref##*/}"_gem_index.gem -i ${species}_mappability_${kmer}.mappability -o ${species}_mappability_${kmer}
${GEM}/wigToBigWig ${species}_mappability_${kmer}.wig ${species}_mappability_${kmer}.sizes ${species}_mappability_${kmer}.bw
${GEM}/bigWigToBedGraph ${species}_mappability_${kmer}.bw ${species}_mappability_${kmer}.bedGraph
## keeping only sites with mappability of 1 (unique) 
${GEM}/bedGraphTobed ${species}_mappability_${kmer}.bedGraph ${species}_mappability_${kmer}.bed 0.99

## merging overlapping intervals
## using the tool download from here: https://github.com/evodify/genotype-files-manipulations/blob/master/combine_overlapping_BEDintervals.py
python ${GEM}/combine_overlapping_BEDintervals.py -i ${species}_mappability_${kmer}.bed -o ${species}_mappability_${kmer}.merged.bed -v 0

#### stats

## calculating the proportion of the genome with mappability ==1
## using the bedmap tool from bedops: https://github.com/bedops/bedops
${bedmap-typical} --echo --bases-uniq --delim '\t' ${ref}.bed ${species}_mappability_${kmer}.merged.bed | \
awk 'BEGIN { genome_length = 0; masked_length = 0; } { genome_length += ($3 - $2); masked_length += $4; } \
END { print (masked_length / genome_length); }' - > MAP1prop.txt

## average mappability over scaffolds
cat ${species}_mappability_${kmer}.bedGraph | awk '{ sum[$1]+=$4; cnt[$1]++} END { print "Scaffold" "\t" "sum" "\t" "cnt" "\t" "avg"; \
for (i in sum) print i "\t" sum[i] "\t" cnt[i] "\t" sum[i]/cnt[i]}' > ${species}_mappability_${kmer}_overScaffolds.txt 

## overlapping regions between the ref and the map1 bed
${bedmap-typical} --echo-overlap-size ${ref}.bed ${species}_mappability_${kmer}.merged.bed | tr -s ';'  '\n' > bedOverlap.txt
${bedmap-typical} --echo-overlap-size ${ref}.bed ${species}_mappability_${kmer}.merged.bed | tr -s ';'  '\n' | awk '{sum +=$1} END{print sum}' > bedOverlapSum.txt

## inverting the "good" bed into a blacklist bed
${bedtools} subtract -a ${ref}.bed -b ${species}_mappability_${kmer}.merged.bed > ${species}_mappability_${kmer}.merged.blacklist.bed

cat muskox_mappability_150.bedGraph | awk '{ if ($4 == 1.0) {print}}' > allGoodregions.bed
${bedtools} intersect -a ${ref}.bed -b allGoodregions.bed -wo > bedintersect_overlap
less bedintersect_overlap | awk '{sum += $8} END{print sum}'
