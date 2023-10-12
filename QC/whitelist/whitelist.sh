#!/bin/bash

bedtools="bedtools2/bin/bedtools"
refBedFile="muskox_pseudohap.gapcloser.v1.lenchecked.fasta.bed"
maskBad="muskox_pseudohap.gapcloser.v1.lenchecked.fasta.RM.bed"
maskGood="muskox_pseudohap.gapcloser.v1.lenchecked.RM.masked.bed"
depthBad="all_remove.bed"
depthGood="all_keep.bed"
mappabilityGood="allGoodregions.bed"
mappabilityBad="muskox_mappability_150.merged.blacklist.bed"
autosomes1M="contigs1MB_autosomes.bed"
inbreedBad="blacklistInbreedSites_unmerged_unsorted.BED"

## sort the reference bed
less $refBedFile | sort -k1,1 -k2,2n > ref.sorted.bed

## subtract bad regions identified by the GEM mappabiility analysis
${bedtools} subtract -a $autosomes1M -b $mappabilityBad > whitelist_part1.bed

## subtract  bad regions with extreme inbreeding coefficient distribution
${bedtools} subtract -a whitelist_part1.bed -b $inbreedBad > whitelist_part2.bed

## subtract  bad regions with extreme sequencing depth distribution
${bedtools} subtract -a whitelist_part2.bed -b $depthBad > whitelist_part3.bed

## subtract bad regions identified by RepeatMasker from the reference bed
## the whitelist_part0.bed should be identical to the masked reference in the ref folder $maskGood
${bedtools} subtract -a whitelist_part3.bed -b $maskBad > whitelist.bed

## intersect the reference bed with the whitelist bed to identified regions of overlap
${bedtools} intersect -a $refBedFile -b whitelist.bed -wo > whitelist_overlap

## sum the regions of overlap to determine final proportion of ref kept
less whitelist_overlap | awk '{ sum += $7 } END { print sum }' > whitelist_sum

awk '{print $1"\t"$2+1"\t"$3}' whitelist.bed > whitelist.regions 
$angsd sites index whitelist.regions 

## subtract bad regions identified by RepeatMasker from the reference bed
## the whitelist_part0.bed should be identical to the masked reference in the ref folder $maskGood
${bedtools} subtract -a $refBedFile -b $maskBad > whitelist_part0.bed

## subtract bad regions identified by the GEM mappabiility analysis
${bedtools} subtract -a whitelist_part0.bed -b $mappabilityBad > whitelist_part1.bed

## subtract non-autosome and weird scaffolds, and sort by scaffold size
cat $autosomes1M | while read line; do grep -w ^"$line" whitelist_part1.bed; done > whitelist_part2.bed

## subtract  bad regions with extreme sequencing depth distribution
${bedtools} subtract -a whitelist_part2.bed -b $depthBad > whitelist_part3.bed

## subtract  bad regions with extreme inbreeding coefficient distribution
${bedtools} subtract -a whitelist_part3.bed -b $inbreedBad > whitelist.bed

## sort the reference bed 
less $refBedFile | sort -k1,1 -k2,2n > ref.sorted.bed

## intersect the reference bed with the whitelist bed to identified regions of overlap
${bedtools} intersect -a ref.sorted.bed -b whitelist.bed -wo > whitelist_overlap

## sum the regions of overlap to determine final proportion of ref kept
less whitelist_overlap | awk '{ sum += $7 } END { print sum }' > whitelist_sum
