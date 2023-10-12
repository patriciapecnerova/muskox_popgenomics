outdir=results
vcf=merged_snps_highcov_10_3_strictSib_rename_withSheep_muskoxsites_diallelic.bcf.gz
bcftools="bcftools"
python="python3.7"
mkdir -p chroms

$bcftools index --stats $vcf | cut -f 1  > chroms.txt
## watch out for underscores in samples names, they will mess up the later script
$bcftools query -l $vcf > samplenames.txt
cat chroms.txt | while read line;
do
    echo "$bcftools query -r ${line} -f '[%GT ]"\\"\n' ${vcf} | $python parse_query.py samplenames.txt  chroms/${line}"
done > parallel.tasks


parallel -j 40 --joblog parallel.log --progress   :::: parallel.tasks

$python merge_2dsfs.py chroms chroms.txt

Rscript reich_fst.R
