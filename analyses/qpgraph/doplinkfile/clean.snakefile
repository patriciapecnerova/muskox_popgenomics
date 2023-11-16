### messy patchy snakemake, copied rules from different pipelines to adapt
## to specific needs in muskox. Do not use unless you know what is it doing and why


BCFTOOLS="/home/genis/software/bcftools/bcftools"
PLINK="/home/genis/software/plink"
REFFINDER="/home/genis/software/refFinder/refFinder"
BGZIP="bgzip"

OUTMAIN=config["outmain"]

rule all:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz"),
        csi = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz.csi"),
        csi2 = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3_nomissing.bcf.gz.csi"),
        withref = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz"),
        csi3 = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz.csi"),
        csi4 = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_nomissing.vcf.gz.csi")



rule keep_highdep:
    input:
        muskox = config["muskox"],
        keep_samples = config["keep_samples"],
        bedkeep = config["bed"]
    output:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz")
    shell: """
    {BCFTOOLS} view --samples-file {input.keep_samples} -T {input.bedkeep}  {input.muskox} | {BCFTOOLS} view --min-ac 1 -Ob -o {output.bcf}
"""


rule index_highdep:
    input:
        bcf=os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz")
    output:
        csi=os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz.csi")
    shell:
        """{BCFTOOLS} index {input.bcf}"""


rule muskox_sites:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz")
    output:
        muskoxsites = temp(os.path.join(OUTMAIN, "bcf", "muskox_sites.txt"))
    shell:
        """{BCFTOOLS} query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.bcf} > {output.muskoxsites}"""


rule merge_bcf:
    input:
        muskox = os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz"),
        muskox_idx = os.path.join(OUTMAIN, "bcf","highdep_snps_10_3.bcf.gz.csi"),
        sheep = config["sheep"],
        muskoxsites = os.path.join(OUTMAIN, "bcf", "muskox_sites.txt")
    output:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz")
    shell:
        """{BCFTOOLS} merge -Ou {input.muskox} {input.sheep} | {BCFTOOLS} view -T {input.muskoxsites} -m2 -M2 -v snps -Ob -o {output.bcf}"""


rule index_bcf_merged:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz")
    output:
        csi = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz.csi")
    shell:
        """{BCFTOOLS} index {input.bcf}"""


rule filter_missing_bcf:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz")
    output:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3_nomissing.bcf.gz")
    shell:
        """{BCFTOOLS} view -e 'F_MISSING > 0' -Ob -o {output.bcf}  {input.bcf}"""


rule index_bcf_nomissing:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3_nomissing.bcf.gz")
    output:
        csi = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3_nomissing.bcf.gz.csi")
    shell:
        """{BCFTOOLS} index {input.bcf}"""



# merged sheep sample
rule bcf_sheep_to_plink:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz"),
        csi = os.path.join(OUTMAIN,"bcf", "highdep_sheep_snps_10_3.bcf.gz.csi")
    output:
        plink = multiext(os.path.join(OUTMAIN, "plink", "merged_sheep_snps_filtered"), ".bed", ".bim", ".fam")
    params:
        prefix = os.path.join(OUTMAIN, "plink", "merged_sheep_snps_filtered"),
    shell:
        """
        {BCFTOOLS} view  -m2 -M2 -v snps {input.bcf} | {BCFTOOLS} annotate -Ov -x ID -I +'%CHROM:%POS:%REF:%ALT' | {PLINK} --vcf /dev/stdin --const-fid --make-bed --allow-extra-chr --out {params.prefix}
"""

        

## convert to plink for admixtools2 and to add ref as sample
rule bcf_to_plink:
    input:
        bcf = os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz"),
        csi = os.path.join(OUTMAIN,"bcf", "highdep_snps_10_3.bcf.gz.csi")
    output:
        plink = multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered"), ".bed", ".bim", ".fam")
    params:
        prefix = os.path.join(OUTMAIN, "plink", "merged_snps_filtered"),
    shell:
        """
        {BCFTOOLS} view  -m2 -M2 -v snps {input.bcf} | {BCFTOOLS} annotate -Ov -x ID -I +'%CHROM:%POS:%REF:%ALT' | {PLINK} --vcf /dev/stdin --const-fid --make-bed --allow-extra-chr --out {params.prefix}
"""


rule bed_to_tped:
    input:
        multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered"), ".bed", ".bim", ".fam")
    output:
        temp(multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered"), ".tped", ".tfam"))
    params:
        prefix = os.path.join(OUTMAIN, "plink", "merged_snps_filtered")
    shell: """
    {PLINK} --bfile {params.prefix} --recode transpose --out {params.prefix} --allow-extra-chr
"""


rule get_ref_tped:
    input:
        bim = os.path.join(OUTMAIN, "plink", "merged_snps_filtered.bim"),
        ref = config["ref"]
    output:
        tped = os.path.join(OUTMAIN, "plink", "ref.tped"),
        tfam = os.path.join(OUTMAIN, "plink", "ref.tfam")
    params:
        os.path.join(OUTMAIN, "plink", "ref")
    shell: """
    cat {input.bim} | {REFFINDER} {input.ref} plink > {output.tped}
    echo "DomesticPig ref 0 0 0 -9" > {output.tfam}
"""


rule add_ref:
    input:
        my_tped = os.path.join(OUTMAIN, "plink", "merged_snps_filtered.tped"),
        my_tfam = os.path.join(OUTMAIN, "plink", "merged_snps_filtered.tfam"),
        ref_tped = os.path.join(OUTMAIN, "plink", "ref.tped"),
        ref_tfam = os.path.join(OUTMAIN, "plink", "ref.tfam"),
    output:
        tped = os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref.tped"),
        tfam = os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref.tfam")
    shell: """
    paste -d" " {input.my_tped} {input.ref_tped} > {output.tped}
    cat {input.my_tfam} {input.ref_tfam} > {output.tfam}
"""


rule tped_to_bed_withref:
    input:
        multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref"),".tped", ".tfam")  
    output:
        multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref"),".bed", ".fam", ".bim")
    params:
        prefix = os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref")
    shell:"""
    {PLINK} -tfile {params.prefix} --make-bed --out {params.prefix} --allow-extra-chr
"""



rule bed_to_vcf_withref:
    input:
        multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref"),".bed", ".fam", ".bim"),
    output:
        os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref.vcf")
    params:
        inprefix = os.path.join(OUTMAIN, "plink", "merged_snps_filtered_withref"),
        outprefix = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref")
    shell:"""
    {PLINK} -bfile {params.inprefix} --recode vcf --out {params.outprefix} --allow-extra-chr
"""



rule fix_ref_alt_vcf:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref.vcf")
    output:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz")
    params:
        ref=config["ref"]
    shell:"""
    {BCFTOOLS} norm -c s -f {params.ref} -N -Oz -o {output.vcf} {input.vcf}
"""


rule index_fixed_ref_alt_vcf:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz")
    output:
        csi = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz.csi")
    shell:
        "{BCFTOOLS} index {input.vcf}"


rule fixed_ref_alt_vcf_nomissing:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz")
    output:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_nomissing.vcf.gz")
    params:
        ref=config["ref"]
    shell:"""
    {BCFTOOLS} view -e 'F_MISSING>0' -Oz -o {output.vcf} {input.vcf}
"""


rule index_fixed_ref_alt_vcf_nomissing:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_nomissing.vcf.gz")
    output:
        csi = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_nomissing.vcf.gz.csi")
    shell:
        "{BCFTOOLS} index {input.vcf}"
