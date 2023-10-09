library(RcppCNPy)

estF <- npyLoad("muskox_whitelist.inbreed.sites.npy")
site <- scan("muskox_whitelist.sites",what="df")
lrt <- npyLoad("muskox_whitelist.lrt.sites.npy")
chr <- sub("_[0123456789]+$","",x=site)
pos <- as.integer(sub(".*_","",x=site))

chr1mb <- scan("contigs1MB_remap_autosomes_headers.txt", what="ds")
k <- chr %in% chr1mb

estF <- estF[k]
lrt <- lrt[k]
chr <- chr[k]
pos <- pos[k]
maxPos <- tapply(pos,chr,max)

meanF <- tapply(estF,chr,mean)
lenF <- tapply(estF,chr,length)

badRegions <- function(ch, minF=-0.95, reg=50000){
    k <- chr == ch
    # get bad sites (meanF below minF and significantly deviating hwe)
    badpos <- pos[k][lrt[k] > 24 & estF[k] < minF]
    
    if(length(badpos)==0) return(data.frame(chr=ch, start=1, end=1))
    # get regions reg (default 50000) bp in both directions of bad sites                             
    badreg <- matrix(c(badpos-reg, badpos+reg), ncol=2)
    badreg[badreg<0] <- 1

    if(nrow(badreg)==1) return(data.frame(chr=ch,start= badreg[1,1], end= badreg[1,2]))
    # collapse overlapping regions
    badreg2 <- c()
    start <- badreg[1,1]
    end <- badreg[1,2]
    for(i in 2:nrow(badreg)){

        if(badreg[i,1]<end){
            end <- badreg[i,2]
        }else{
            badreg2 <- c(badreg2, start, end)
            start <- badreg[i,1]
            end <- badreg[i,2]
        }
    }

    badreg2 <- c(badreg2,start,end)
    
    badreg <- t(matrix(badreg2,nrow=2))
    #out <- cbind(rep(ch, nrow(badreg)), badreg)
    out <- data.frame(chr=rep(ch,nrow(badreg)), start = badreg[,1], end = badreg[,2])
    return(out) 
}



badbedl <- lapply(chr1mb, badRegions,minF=-0.95, reg=50000)

badbed <- do.call('rbind', badbedl)

allbed <- read.table("contigs1MB_remap_autosomes.bed")


names(allbed) <- names(badbed)

lenghtbad <- tapply(badbed$end-badbed$start, badbed$chr, sum)

Nbadreg <- tapply(badbed$chr, badbed$chr, length)* (lenghtbad>0)

summarydf <- data.frame(chr=chr1mb, length=allbed$end,meanF=meanF[as.factor(chr1mb)],
                        proportionBad=lenghtbad/allbed$end,
                        Nbadregions=Nbadreg,
                        keep= !(lenghtbad/allbed$end > 0.2 | meanF < -0.2))
write.table(summarydf, "inbreedSummary.tsv",
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


w <- order(meanF)
pdf("inbreedfilts1mb.pdf")
for(i in 1:length(w)){
     plot(pos[chr==names(meanF[as.factor(chr1mb)])[w[i]]],estF[chr==names(meanF[as.factor(chr1mb)])[w[i]]],pch=4,col="goldenrod",lwd=2,main=names(meanF[as.factor(chr1mb)])[w[i]],ylab="F",xlab="position")
         p<-pos[chr==names(meanF[as.factor(chr1mb)])[w[i]]]
         F<-estF[chr==names(meanF[as.factor(chr1mb)])[w[i]]]
         win <- 200
         fun<-function(x,w)
            ( cumsum(as.numeric(x))[-c(1:w)]-rev(rev(cumsum(as.numeric(x)))[-c(1:w)])) / w
    lines(fun(p,win),fun(F,win))
     abline(v=badbed[badbed$chr==chr1mb[w[i]],]$start, col="darkred",lwd=3,lty=2)
     abline(v=badbed[badbed$chr==chr1mb[w[i]],]$end, col="darkred",lwd=3,lty=2)
     }
dev.off()


badbed <- badbed[badbed$end > 1,]

finalbadbed <- rbind(badbed, data.frame(chr=summarydf$chr[!summarydf$keep],
                                start=rep(1, sum(!summarydf$keep)),
                                end=summarydf$length[!summarydf$keep]))

write.table(finalbadbed, "blacklistInbreedSites_unmerged_unsorted.BED", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# this create blacklist, which then is merged and used to create whitelist with bedtools
