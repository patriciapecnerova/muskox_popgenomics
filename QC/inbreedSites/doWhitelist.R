## create BED

blacklist <- read.table("blacklistInbreedSites_unmerged_unsorted.REGIONS")
blacklist2 <- blacklist
blacklist2[,2:3] <- blacklist[,2:3] - 1
write.table(blacklist2, "blacklistInbreedSites_unmerged_unsorted.BED", col.names=F, row.names=F, quote=F, sep="\t")

## turn the BED into angsd
bedall <- read.table("contigs1MB_remap_autosomes.bed")
bedall2 <- bedall
bedall2[,2:3] <- bedall[,2:3] + 1
write.table(bedall2, "contigs1MB_remap_autosomes.angsd", col.names=F, row.names=F, quote=F, sep="\t")
