
BCFTOOLS="bcftools"
SETGT="bcftools/plugins/setGT.so"
MIN_MQ=30
MIN_BQ=30
OUTMAIN=config["outmain"] # "/emc/kristian/leopard/vcfs/results"
REF=config["ref"]

configfile: "config_vcf_muskoxmap_Sib.yaml"

## this is a filter used together with total depth and allele depth
BED = config["bed"]

print(OUTMAIN)
print(REF)
print(BED)

#### defaults
depths = ["10", "5"]
alleledepths = ["2", "3"]

wildcard_constraints:
    sample = "|".join(config["samples"].keys()),
    chrom="|".join(config["chroms"]),
    dp = "|".join(depths),
    ad = "|".join(alleledepths),


rule all:
    input:
        expand(os.path.join(OUTMAIN, "vcf", "{sample}.bcf.gz"), sample=config["samples"].keys()),
        os.path.join(OUTMAIN, "vcf_highcov", "merged_snps_highcov.bcf.gz"),
        expand(os.path.join(OUTMAIN, "vcf_highcov_filtered", "merged_snps_highcov_{dp}_{ad}.bcf.gz"),
               dp=depths,
               ad=alleledepths),
        expand(os.path.join(OUTMAIN, "vcf_filtered", "{sample}_{dp}_{ad}.bcf.gz"),
               sample=config["samples"].keys(),
               dp=depths,
               ad=alleledepths),


    shell:
        "ln -s -f {OUTMAIN} ."


rule gen_bcftools_genome_wide_indi:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp(os.path.join(OUTMAIN, "vcf", "{sample}_{chrom}.bcf.gz"))
    threads: 1
    shell:
        """{BCFTOOLS} mpileup -r {wildcards.chrom} -G readgroup.txt -B -Q {MIN_BQ}  -q {MIN_MQ} --threads {threads} -O u --fasta-ref {REF} --per-sample-mF -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR {input}  | {BCFTOOLS} call -Ob -o {output} --threads {threads} --multiallelic-caller
{BCFTOOLS} index {output}
"""


rule concat_chroms_indi:
    input:
        expand(os.path.join(OUTMAIN, "vcf", "{{sample}}_{chrom}.bcf.gz"), chrom = config["chroms"])
    output:
        vcf=os.path.join(OUTMAIN, "vcf", "{sample}.bcf.gz"),
        idx=os.path.join(OUTMAIN, "vcf", "{sample}.bcf.gz.csi")
    shell: """
    {BCFTOOLS} concat --naive -Ob -o {output.vcf} {input}
    {BCFTOOLS} index {output.vcf}
    """

rule filter_indi:
    input:
        vcf=os.path.join(OUTMAIN, "vcf", "{sample}.bcf.gz"),
    output:
        vcf=os.path.join(OUTMAIN, "vcf_filtered", "{sample}_{dp}_{ad}.bcf.gz"),
        idx=os.path.join(OUTMAIN, "vcf_filtered", "{sample}_{dp}_{ad}.bcf.gz.csi"),
    threads: 3
    shell: """
    {BCFTOOLS} +{SETGT} --threads {threads} -T {BED} -Ob -o {output.vcf} {input.vcf} -- -t q -n . -i  '(GT="hom" & FMT/DP<{wildcards.dp}) | (GT="het" & (FMT/AD[:0]<{wildcards.ad} | FMT/AD[:1]<{wildcards.ad} | FMT/DP<{wildcards.dp}))'
    {BCFTOOLS} index {output.vcf}
    """


rule merge_indi:
    input:
        expand(os.path.join(OUTMAIN, "vcf", "{sample}.bcf.gz"), sample=config["samples"].keys())
    output:
        vcf=os.path.join(OUTMAIN, "vcf_highcov", "merged_snps_highcov.bcf.gz"),
        # csi=os.path.join(OUTMAIN, "vcf_highcov", "merged_snps_highcov.bcf.gz.csi")
    shell: """
    {BCFTOOLS} merge --force-samples --threads 24 -Ou {input} | {BCFTOOLS} view --threads 24 -v snps -Ob -o {output.vcf}
    # {BCFTOOLS} index {output.vcf}
    """

### IMPORTANT & and | means within same sample. && and || checks all the samples.

rule filter_merged:
    input:
        vcf=os.path.join(OUTMAIN, "vcf_highcov", "merged_snps_highcov.bcf.gz"),
    output:
        vcf=os.path.join(OUTMAIN, "vcf_highcov_filtered", "merged_snps_highcov_{dp}_{ad}.bcf.gz"),
        idx=os.path.join(OUTMAIN, "vcf_highcov_filtered", "merged_snps_highcov_{dp}_{ad}.bcf.gz.csi"),
    threads: 8
    shell: """
    {BCFTOOLS} +{SETGT} --threads {threads} -T {BED} -Ob -o {output.vcf} {input.vcf} -- -t q -n . -i  '(GT="hom" & FMT/DP<{wildcards.dp}) | (GT="het" & (FMT/AD[:0]<{wildcards.ad} | FMT/AD[:1]<{wildcards.ad} | FMT/DP<{wildcards.dp}))'
    {BCFTOOLS} index {output.vcf}
    """
