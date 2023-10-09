BEDTOOLS="bedtools2/bin/bedtools"

bedall="contigs1MB_remap_autosomes.bed"

blacklist="blacklistInbreedSites"
#whitelist=whitelistInbreedSites

$BEDTOOLS sort -i ${blacklist}_unmerged_unsorted.BED > ${blacklist}_unmerged.BED
$BEDTOOLS merge -i ${blacklist}_unmerged.BED > ${blacklist}.BED

$BEDTOOLS subtract -a $bedall -b $blacklist.BED > $whitelist.BED



### R Check % sequence in blacklist / whitelist

whitelist <- read.table("whitelistInbreedSites_unmerged_unsorted.BED")
blacklist <- read.table("blacklistInbreedSites_unmerged_unsorted.BED")

alllist <- read.table("contigs1MB_remap_autosomes.bed")

lenlist <- function(x) sum(as.numeric(x$V3 - x$V2))

lenlist(whitelist) / lenlist(alllist) # 0.9480168
lenlist(blacklist) / lenlist(alllist) # 0.05279217


lenlist(whitelist) / lenlist(alllist) + lenlist(blacklist) / lenlist(alllist) # 1.000809 above 1 bc some regions in blacklist span longer than chromosome length


#### Make bed files 0 indexed instead of 1 indexed (by subracting 1 from each value)
#write.table(whitelist2, "whitelistInbreedSites.BED", col.names=F, row.names=F, quote=F, sep="\t")

#write.table(blacklist2, "blacklistInbreedSites.BED", col.names=F, row.names=F, quote=F, sep="\t")


whitelist <- read.table("whitelistInbreedSites.REGIONS")
blacklist <- read.table("blacklistInbreedSites.REGIONS")
alllist <- read.table("GCF_002742125.1_Oar_rambouillet_v1.0_genomic_whitelist_autosomes.ANGSD")

lenlist <- function(x) sum(as.numeric(x$V3 - x$V2))

lenlist(whitelist) / lenlist(alllist) # 0.9480168
lenlist(blacklist) / lenlist(alllist) # 0.05279217


lenlist(whitelist) / lenlist(alllist) + lenlist(blacklist) / lenlist(alllist) # 1.000809 above 1 bc some regions in blacklist span longer than chromosome length


#### Make bed files 0 indexed instead of 1 indexed (by subracting 1 from each value)
whitelist2 <- whitelist
whitelist2[,2:3] <- whitelist[,2:3] - 1
write.table(whitelist2, "whitelistInbreedSites.BED", col.names=F, row.names=F, quote=F, sep="\t")

blacklist2 <- blacklist
blacklist2[,2:3] <- blacklist[,2:3] - 1
write.table(blacklist2, "blacklistInbreedSites.BED", col.names=F, row.names=F, quote=F, sep="\t")
