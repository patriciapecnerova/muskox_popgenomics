#!/bin/bash

evalAdmix="evalAdmix"
ref="muskox_pseudohap.gapcloser.v1.lenchecked.fasta"
beagle="muskox_muskoxmap_whitelist.beagle.gz"
output="results"
admixResultFolder="results"

$evalAdmix/evalAdmix -beagle $beagle -fname "$admixResultFolder"/5/ngsAdmix_muskox_remap2.5.4.fopt.gz -qname "$file".qopt -P 2 -o $output evaladmix_muskoxmap_remap2_convergence_K5

$evalAdmix/evalAdmix -beagle $beagle -fname "$admixResultFolder"/6/ngsAdmix_muskox_remap2.6.6.fopt.gz -qname "$file".qopt -P 2 -o $output evaladmix_muskoxmap_remap2_convergence_K6

$evalAdmix/evalAdmix -beagle $beagle -fname "$admixResultFolder"/7/ngsAdmix_muskox_remap2.7.2.fopt.gz -qname "$file".qopt -P 2 -o $output evaladmix_muskoxmap_remap2_convergence_K7

$evalAdmix/evalAdmix -beagle $beagle -fname "$admixResultFolder"/8/ngsAdmix_muskox_remap2.8.8.fopt.gz -qname "$file".qopt -P 2 -o $output evaladmix_muskoxmap_remap2_convergence_K8

$evalAdmix/evalAdmix -beagle $beagle -fname "$admixResultFolder"/9/ngsAdmix_muskox_remap2.9.28.fopt.gz -qname "$file".qopt -P 2 -o $output evaladmix_muskoxmap_remap2_convergence_K9
