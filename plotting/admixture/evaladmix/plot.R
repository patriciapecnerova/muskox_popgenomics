source("/home/genis/github/evalAdmix/visFuns.R")

# read population labels and estimated admixture proportions
popinfo<-read.table("het_remap2_wSib.txt",as.is=T)

region <- popinfo$V1
inds <- popinfo$V2
loc <- paste(region, gsub("EL", "SIB", gsub("MO", "ZA", substr(inds,1,2))), sep="_")


popord <- c("SIB", "CaMW", "CaME", "CaIS", "CaIN", "GrN", "GrEN", "GrES")
ord <- orderInds(pop=region, popord=popord)
locord <- unique(loc[ord])

##### K6
q <- as.matrix(read.table("ngsAdmix_muskox_remap2_wSib.6.5.qopt"))

ord <- orderInds(pop=loc, q=q, popord=locord)

#ord <- orderInds(q, pop = region, popord = popord)
refinds <- c('CM-1', 'HD-5', 'MO32', 'AH-1', 'KU-16', 'CH-4')
refidx <- sapply(refinds, function(x) which(inds==x))
q6 <- q[,orderK(q, refinds = refidx[1:6])]

muskoxRegCol <- c("CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853",
                  "GrEN" = "#D02A7B", "GrES" = "#F17A30")
colorord <- c(3,6,5,4,1,2)
colorpal <- muskoxRegCol[colorord]


barplot(t(q6)[,ord],col=colorpal,space=0,border =NA, 
        cex.names=1, cex.axis = 2.3, cex.lab = 2.5,
        ylab="Admixture proportions for K=6", names.arg = inds[ord], las=2)
text(tapply(1:length(region),region[ord],mean),-0.32,unique(region[ord]),xpd=T, cex=2)
abline(v=cumsum(sapply(unique(region[ord]),function(x){sum(region[ord]==x)})),col=1,lwd=1.4)
mtext("Individuals", side=1, line=3, cex=2.3, padj = 4)


loclab <- unlist(sapply(strsplit(loc, "_"), function(x) x[2]))

corres <- as.matrix(read.table("evaladmix_muskoxmap_remap2_wSib_convergence_K6"))

for(k in 5:8){
    f <- paste0("evaladmix_muskoxmap_remap2_wSib_convergence_K", k)
    outpng <- paste0("muskox_corresK", k, ".png")
    corres <- as.matrix(read.table(f))
    bitmap(outpng, w=6.5, h=6,res=300)
    plotCorRes(cor_mat=corres, pop=loclab, ord=ord, superpop=region, max_z=0.1, cex.lab=0.9, cex.lab.2=1,rotatelabpop=90, title=paste0("Evaluation of muskox admixture proportions with K=",k), jitterpoplab=TRUE)
    dev.off()
}
        
