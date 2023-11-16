library(admixtools)

outtable <- "/home/genis/other/muskox/dstats_indiv/results/muskoxindiv_dstat_sibp3.tsv"
plFile <- "/home/genis/other/muskox/dstats_indiv/indata/merged_snps_filtered_withref"

fam <- read.table(paste0(plFile, ".fam"))


# p1 and p2 will be all samples except siberian and sheep
p1 <- fam$V2[2:108]
p2 <- fam$V2[2:108]

p3 <- "EL012"
p4 <- "SheepRef"

res <- qpdstat(plFile, pop1=p1, pop2=p2, pop3=p3, pop4=p4,
               f4Mode=FALSE, blgsize=5e6,unique_only=TRUE,
               maxmiss=0, allsnps=TRUE, auto_only = FALSE)


write.table(res, outtable, col.names=T, row.names=F, quote=F, sep="\t")

cat("FINISHED AND WROTE DSTATS TO", outtable,"\n")

