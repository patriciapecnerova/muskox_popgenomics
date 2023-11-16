# Rscript extractf2s.R inplink popinfo outdir blocklen maxmis

library(admixtools)

args <- commandArgs(trailingOnly=T)

inplink <- args[1]
popinfo <- args[2]
f2dir <- args[3]
blocklen <- args[4]
maxmiss <- args[5]

info <- read.table(popinfo,h=F)

pops <- info$V1
inds <- info$V2

extract_f2(
    pref=inplink,
    outdir=f2dir,
    pops=pops,
    inds=inds,
    maxmiss=maxmiss,
    format="plink",
    blgsize=blocklen,
    fst=TRUE
)
