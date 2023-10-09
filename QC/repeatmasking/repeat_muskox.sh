#!/bin/bash

RepeatModeler="RepeatModeler-2.0.1"
RepeatMasker="RepeatMasker"
REF="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
out="results"
references="references"

$RepeatModeler/BuildDatabase -name muskox $REF -engine rmblast

$RepeatModeler/RepeatModeler -engine rmblast -pa 20 -database muskox > muskox_RepeatModeler.log

perl $RepeatMasker/RepeatMasker -pa 8 -a -gccalc -dir ./ -lib $out/consensi.fa.classified $REF

python3.5 $RepeatMasker/util/RM2Bed.py muskox_pseudohap.gapcloser.v1.lenchecked.fasta.out

bedtools2/bin/bedtools subtract -a $references/muskox_pseudohap.gapcloser.v1.lenchecked.fasta.bed -b muskox_pseudohap.gapcloser.v1.lenchecked.fasta.RM.bed > muskox_pseudohap.gapcloser.v1.lenchecked.RM.masked.bed
