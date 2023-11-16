library(admixtools)

f2dir <- "/home/genis/other/muskox/qpgraph2_sheepmap/f2_statistics"
plFile <- "/home/genis/other/muskox/qpgraph2_sheepmap/plink/merged_snps_filtered_withref"

fam <- read.table(paste0(plFile,".fam"), h=F, stringsAsFactors=F)
#pops <- unique(fam$V1)

popfile <- "/home/genis/other/muskox/qpgraph2_sheepmap/pop.info"
popinfo <- read.table(popfile, stringsAsFactors=F, h=F)

pops <- popinfo$V2
inds <- fam$V1

## first, run one time, extract f2 values between pops and save in disk:
extract_f2(pref=plFile, outdir=f2dir, pops=pops, inds=inds, maxmiss=0.1, format="plink",
           blgsize=5e6, fst=TRUE, auto_only=FALSE)


#f2s <- read_f2(f2dir)

