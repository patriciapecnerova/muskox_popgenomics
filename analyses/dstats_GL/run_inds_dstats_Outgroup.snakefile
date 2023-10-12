import pandas as pd
import re

ANGSD="angsd"
JACKKNIFE="angsd/R/jackKnife.R"

OUTMAIN="results_all"

configfile:"config_wSib_wSheep.yaml"

infofile=config['info']

info = pd.read_table(infofile)



wildcard_constraints:
    group = "|".join(config["groups"].keys())



rule all:
    input:
        expand(os.path.join(OUTMAIN, "dstats", "{group}_dstats.txt"), group=config["groups"].keys())


rule do_bamfiles:
    input:
        config['info']
    output:
        bamlists = expand(os.path.join(OUTMAIN, "bamlists", "{group}.bamlist"), group=config['groups'].keys()),
        ids =  expand(os.path.join(OUTMAIN, "bamlists", "{group}.ids"), group=config['groups'].keys())
    run:
        for group in config['groups'].keys():
            out_bam = os.path.join(OUTMAIN, "bamlists", "{}.bamlist").format(group)
            out_ids = os.path.join(OUTMAIN, "bamlists", "{}.ids").format(group)
            info.Bam.loc[info.Region.apply(lambda x: x in config['groups'][group])].to_csv(out_bam, index=False, header=False)
            info.Sample.loc[info.Region.apply(lambda x: x in config['groups'][group])].to_csv(out_ids, index=False, header=False)


rule do_abbababa_per_chrom:
    input:
        bamlist = os.path.join(OUTMAIN, "bamlists", "{group}.bamlist")
    output:
        temp(os.path.join(OUTMAIN, "abbababa", "{group}_{chrr}.abbababa"))
    params:
        outprefix = os.path.join(OUTMAIN, "abbababa", "{group}_{chrr}"),
        r = "{chrr}",
        regions = config['bed'],
        MINQ = 30,
        MINMAPQ = 30,
        blocksize = 5000000,
        anc = config['anc']
    log: os.path.join(OUTMAIN, "abbababa", "{group}_{chrr}.arg")
    shell: "{ANGSD} -doAbbababa 1 -doCounts 1 -P 10 -GL 2 -out {params.outprefix}  -bam {input.bamlist} -minmapQ {params.MINMAPQ} -minq {params.MINQ} -r {params.r} -sites {params.regions} -blockSize {params.blocksize} -anc {params.anc}"




rule concat_abbababa:
    input:
        expand(os.path.join(OUTMAIN, "abbababa", "{{group}}_{chrr}.abbababa"), chrr=config["chroms"])
    output:
        f = os.path.join(OUTMAIN, "abbababa", "{group}.abbababa")
    run:
        f1 = input[0]
        shell("head -1 {f1} > {output.f}")
        for f in input:
            shell("cat {f} | sed 1d >> {output.f}")



rule do_jackknife:
    input:
        abbababa = os.path.join(OUTMAIN, "abbababa", "{group}.abbababa"),
        ids = os.path.join(OUTMAIN, "bamlists", "{group}.ids")
    output:
        os.path.join(OUTMAIN, "dstats", "{group}_dstats.txt")
    params:
        outprefix = os.path.join(OUTMAIN, "dstats", "{group}_dstats")
    shell: "Rscript {JACKKNIFE} file={input.abbababa} indNames={input.ids} outfile={params.outprefix}"
