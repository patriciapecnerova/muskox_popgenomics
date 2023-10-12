SAMTOOLS="samtools"
AWK_STUFF = "a.awk"

MY_SITES = "whitelist.bed"

BAM_IN_DIR="mapping2/results/MQ30bam"

configfile: "config_psmc_muskox.yaml"

#SAMPLES = glob_wildcards(os.path.join(BAM_IN_DIR, "{s}.MQ30.bam")).s
#sample = config["samples"]

wildcard_constraints:
    sample = "|".join(config["samples"].keys()),

OUTMAIN="results"

rule all:
     input:
        os.path.join(OUTMAIN, 'all_depths.txt')


rule do_depth_samtools:
     input:
        bam= lambda wildcards: config["samples"][wildcards.sample]
#        bam=expand(os.path.join(BAM_IN_DIR, "{s}.MQ30.bam"), s=SAMPLES)
     output:
        os.path.join(OUTMAIN, "samples", "{sample}.depth")
     params:
#        sample = {s}
     shell: "{SAMTOOLS} depth -a -Q 30 -b {MY_SITES} {input.bam} | awk -f {AWK_STUFF} > {output}"


rule collect_res:
     input:
        expand(os.path.join(OUTMAIN, "samples", "{sample}.depth"), sample=config["samples"].keys())
     output:
        os.path.join(OUTMAIN, "all_depths.txt")
     shell: "cat {input} > {output}"
