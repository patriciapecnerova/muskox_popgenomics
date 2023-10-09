# config needs:
#    samples: dic with pairs of sample:pathtobamfile
#    depths: dic with depth thresholds specific to each sample sample:[mindepth, maxdepth]
#    beds: dic with filtername:path_to_bed_with_filters
#    allele_support: list of allele support values to call heterozygotes
#    plot_params: dic with g:generation_time mu:mutation_rate
#    outmain: path to main output directory

import os

# https://github.com/lh3/psmc

VCFUTILS="bcftools/misc/vcfutils.pl"
SAMTOOLS="samtools"
BCFTOOLS="bcftools"
SETGT="bcftools/plugins/setGT.so"
MIN_MQ=30
MIN_BQ=30
OUTMAIN=config["outmain"]
REF=config["ref"]
BEDS = config["beds"]
allele_support = config["allele_support"]
configfile: "config_het_from_gc.yaml"

## this is a filter used together with total depth and allele depth
print(OUTMAIN)
print(REF)

wildcard_constraints:
    t = "|".join(allele_support),
    sample = "|".join(config["samples"].keys()),
    bed = "|".join(BEDS.keys()),


rule all:
    input:
        expand(os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),

# rules to estimate genome_wide heterozygoisties with bcftools, using same sites as in psmc

rule get_het_mindepthx:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f rm_indels.awk | {BCFTOOLS} stats -s - > {output}
        """

rule estimate_het:
    input:
       os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats")
    output:
       os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het")
    shell: "grep '^PSC' {input} | awk '{{print $6/($4+$5+$6)}}' > {output}"
