
source("/home/genis/github/evalAdmix/visFuns.R")

popinfof <- "/home/genis/other/muskox/plotadmixture/adm_final/het_remap2_wSib.txt"
indir <- "/home/genis/other/muskox/plotadmixture/adm_final"

muskoxcolors <- c("CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853", "GrEN" = "#D02A7B", "GrES" = "#F17A30", "#d0a304", "#919FA9")


muskoxcolors <- c("KU-2" = "#8444DE",
                  "BL-13" = "#618DFF",
                  "CM-1" = "#41A99C",
                  "DI-3" = "#368853",
                  "MO32" = "#D02A7B",
                  "HD-3" = "#F17A30",
                  "CH-3" = "#d0a304",
                  "DI-5" = "#919FA9")


q <- as.matrix(read.table("/home/genis/other/muskox/plotadmixture/adm_final/ngsAdmix_muskox_remap2_wSib.8.3.qopt"))
popinfo <- read.table(popinfof, as.is=T)
region <- popinfo$V1
inds <- popinfo$V2
loc <- paste(region, gsub("EL", "SIB", gsub("MO", "ZA", substr(inds,1,2))), sep="_")
popord <- c("SIB", "CaMW", "CaME", "CaIS", "CaIN", "GrN", "GrEN", "GrES")
ord <- orderInds(pop=region, popord=popord)
locord <- unique(loc[ord])

loclab <- unlist(sapply(strsplit(loc, "_"), function(x) x[2]))

ord <- orderInds(pop=loc, q=q, popord=locord)



inds <- popinfo$V2
refinds <- c('CM-1', 'KU-2','BL-13','DI-3', 'MO32','HD-3','CH-3','DI-5') # CONTINUE DI3
refidx <- sapply(refinds, function(x) which(inds==x))


bitmap("/home/genis/other/muskox/plotadmixture/muskoxadmixk5to8.png", width=8, height=5, res=300)

par(mfrow=c(4,1), mar = c(1, 4, 1, 5), oma=c(5,0,0,0))

for(k in 5:8){
    
    qfile <- list.files("/home/genis/other/muskox/plotadmixture/adm_final", pattern= paste0(k,".[1-9].qopt"), full.names=T)
    q <- as.matrix(read.table(qfile))
    #plotAdmix(q[,orderK(q=q, refpops=kpopord[1:k], pop=info$Locality)],ord=ord, main="",colorpal=impala_colors[kpopord[1:k]], cex.laby=0.9)
    kord <- orderK(q, refinds = refidx[1:k])

    barplot(t(q)[kord,ord], col=muskoxcolors[refinds[1:k]], space=0, border=NA, cex.axis=1.2,cex.lab=1,
            ylab="", xlab="", main="", cex.main=1.5,xpd=NA, yaxt="n")
            
    axis(2, line=-1, cex.axis=1.5)
    title(ylab="Admixture proportion", line=2, xpd=NA, cex.lab=1.75)
#    abline(v=1:nrow(q), col="white", lwd=0.01)


    abline(v=cumsum(sapply(unique(loc[ord]),function(x){sum(loc[ord]==x)}))[-length(unique(loc))],col=1,lwd=1)
    abline(v=cumsum(sapply(unique(region[ord]),function(x){sum(region[ord]==x)}))[-length(unique(region))],col=1,lwd=2)
    text(x=112, y=0.5, labels=paste("K =", k), xpd=NA,cex=2)
    
}

text(x=1:nrow(q)-0.5, y=-0.17,labels= inds[ord],srt=90, xpd=NA,cex=0.9)

text(sort(tapply(1:length(region),region[ord],mean))-1,-0.4,
     unique(region[ord]), xpd=NA, cex=1.5,
     font=2)

dev.off()




k <- 6
bitmap("/home/genis/other/muskox/plotadmixture/muskoxadmixk6.png", width=8, height=3, res=300)

#par( oma=c(2,0,0,0))

qfile <- list.files("/home/genis/other/muskox/plotadmixture/adm_final", pattern= paste0(k,".[1-9].qopt"), full.names=T)
q <- as.matrix(read.table(qfile))
#plotAdmix(q[,orderK(q=q, refpops=kpopord[1:k], pop=info$Locality)],ord=ord, main="",colorpal=impala_colors[kpopord[1:k]], cex.laby=0.9)
kord <- orderK(q, refinds = refidx[1:k])

barplot(t(q)[kord,ord], col=muskoxcolors[refinds[1:k]], space=0, border=NA, cex.axis=1.2,cex.lab=1,
        ylab="", xlab="", main="", cex.main=1.5,xpd=NA, yaxt="n")
axis(2, line=-1, cex.axis=1.5)
title(ylab="Admixture proportion", line=2, xpd=NA, cex.lab=1.75)
#abline(v=1:nrow(q), col="white", lwd=0.01)


abline(v=cumsum(sapply(unique(loc[ord]),function(x){sum(loc[ord]==x)}))[-length(unique(loc))],col=1,lwd=1)
abline(v=cumsum(sapply(unique(region[ord]),function(x){sum(region[ord]==x)}))[-length(unique(region))],col=1,lwd=2)
text(x=112, y=0.5, labels=paste("K =", k), xpd=NA,cex=2)



text(x=1:nrow(q)-0.5, y=-0.075,labels= inds[ord],srt=90, xpd=NA,cex=0.9)

text(sort(tapply(1:length(region),region[ord],mean))-1,-0.2,
     unique(region[ord]), xpd=NA, cex=1.5,
     font=2)

dev.off()
