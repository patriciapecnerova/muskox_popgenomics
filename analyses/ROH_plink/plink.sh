#!/bin/bash

bfile="merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic.bcf.gz"
vfile="merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic.vcf.gz"

bcftools view -Ov $bfile -o $vfile

## 1mb
plink -vcf $vfile --homozyg --out merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic --aec --homozyg-kb 1000 --homozyg-window-missing 40

## 0.5mb
plink -vcf $vfile --homozyg --out merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic --aec --homozyg-kb 500 --homozyg-window-missing 40

## 2mb
plink -vcf $vfile --homozyg --out merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic --aec --homozyg-kb 2000 --homozyg-window-missing 40

## 5mb
plink -vcf $vfile --homozyg --out merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic --aec --homozyg-kb 5000 --homozyg-window-missing 40

## 10mb
plink -vcf $vfile --homozyg --out merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic --aec --homozyg-kb 10000 --homozyg-window-missing 40
