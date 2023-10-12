library(data.table)

library(grDevices)
library(scales)


# Feed it your dstats and information 
df <- read.table("All.ids", h = T, sep = "\t")

###### CHOOSE LOCATION OR REGION!!!
locIdx <- df$Location
locIdx <-df$Region
######################


names(locIdx) <- df$Sample

#admixed <- c("7946", "5180", "5520")

dstatdf <- fread("All_dstats.txt",
data.table = F, stringsAsFactors = F)[,1:9]


# Function to get D-stats for plotting
getDstat3pop <- function(dstatdf, poptree){
  # extracts all individual dststs that meet poptree, where poptree is a vector c(H1,H2,H3)
  # returns also those that meet c(H2,H1,H3) with Dstat and Z score multiplied by -1 (IS THIS A CORRECT THING TO DO??)
  # return value is a dataframe with same columns as input plus poptree column specifying the population tree
  
#  if(!keepAdmix) dstatdf <- dstatdf[!(dstatdf$H1%in%admixed|dstatdf$H2%in%admixed|dstatdf$H3%in%admixed),]
  
  indTrees <- matrix(as.character(as.matrix(dstatdf[,1:3])), ncol=3)
  
  keep1 <- apply(matrix(locIdx[indTrees], ncol=3),  1, function(x) all(x == poptree))
  keep2 <- apply(matrix(locIdx[indTrees], ncol=3),  1, function(x) all(x == poptree[c(2,1,3)]))
  
  dstat1 <- dstatdf[keep1,]
  dstat2 <- dstatdf[keep2,]
  
  dstat2[,c("Dstat","jackEst", "Z")] <-  - dstat2[,c("Dstat","jackEst", "Z")]
  
  dstat <- rbind(dstat1,dstat2)
  
  dstat$poptree <- paste(poptree, collapse = ".")
  
  return(dstat)
}



# get all combinations with Sib as outgroup 

regions<-levels(factor(df$Region))
regions<-regions[-c(1,2,11,16,18,19)] # Removes Sib, Sheep and GrES
pairs <- t(combn(regions, m=2))
treesSibOut <- cbind(pairs, rep("Sib", nrow(pairs)))

dstatsSibOut <- do.call(rbind,apply(treesSibOut, 1, getDstat3pop, dstatdf=dstatdf))
dstatsSibOut$poptree <- as.factor(dstatsSibOut$poptree, levels = treesSibOut)


#regions<-levels(factor(df$Location))
#regions<-regions[-c(6,8:9)] # Removes Sib, Sheep and GrES
#pairs <- t(combn(regions, m=2))
#treesSibOut <- cbind(pairs, rep("Sib", nrow(pairs)))

#dstatsSibOut <- do.call(rbind,apply(treesSibOut, 1, getDstat3pop, dstatdf=dstatdf))
#dstatsSibOut$poptree <- as.factor(dstatsSibOut$poptree, levels = treesSibOut)


png("testBoxplotDstatsSibH3_noGrES.png",width=800, height = 1600)
par(mar=c(5,14,4,4))
boxplot(dstatsSibOut$Z ~ dstatsSibOut$poptree,
        boxwex=0.5,
        main="All individuals Dstat combinations setting Sib as outroup",
        horizontal=T,las=2,ylab="",ylim=c(-10,10),xlab="Z score",cex.axis=1.2,cex.lab=1.5)
abline(v=c(-3,3),lty=2)
dev.off()


# Produce single plot with 5 panels: CHOOSE/ CREATE YOUR GROUPING MATRIX. GROUP MUST HAVE MORE THAN 1 INDIVIDUAL

## Groups by Locality
group_KU <- matrix(c(rep("KU",2), "TH",
                     rep("KU",2), "BL",
                     rep("KU",2), "CH",
                     rep("KU",2), "FR",
                     rep("KU",2), "AI",
                     rep("KU",2), "CM",
                     rep("KU",2), "BI",
                     rep("KU",2), "CI",
                     rep("KU",2), "DI",
                     rep("KU",2), "GF",
                     rep("KU",2), "SF",
                     rep("KU",2), "AH",
                     rep("KU",2), "EU",
                     rep("KU",2), "EI",
                     rep("KU",2), "RY",
                     rep("KU",2), "ZA",
                     rep("KU",2), "JA",
                     rep("KU",2), "HD"),
                   ncol=3, byrow=T)

group_TH <- matrix(c(rep("TH",2), "KU",
                     rep("TH",2), "BL",
                     rep("KU",2), "CH",
                     rep("TH",2), "FR",
                     rep("TH",2), "AI",
                     rep("TH",2), "CM",
                     rep("TH",2), "BI",
                     rep("TH",2), "CI",
                     rep("TH",2), "DI",
                     rep("TH",2), "GF",
                     rep("TH",2), "SF",
                     rep("TH",2), "AH",
                     rep("TH",2), "EU",
                     rep("TH",2), "EI",
                     rep("TH",2), "RY",
                     rep("TH",2), "ZA",
                     rep("TH",2), "JA",
                     rep("TH",2), "HD"),
                   ncol=3, byrow=T)

group_BL <- matrix(c(rep("BL",2), "KU",
                     rep("BL",2), "TH",
                     rep("KU",2), "CH",
                     rep("BL",2), "FR",
                     rep("BL",2), "AI",
                     rep("BL",2), "CM",
                     rep("BL",2), "BI",
                     rep("BL",2), "CI",
                     rep("BL",2), "DI",
                     rep("BL",2), "GF",
                     rep("BL",2), "SF",
                     rep("BL",2), "AH",
                     rep("BL",2), "EU",
                     rep("BL",2), "EI",
                     rep("BL",2), "RY",
                     rep("BL",2), "ZA",
                     rep("BL",2), "JA",
                     rep("BL",2), "HD"),
                   ncol=3, byrow=T)

group_CH <- matrix(c(rep("CH",2), "KU",
                     rep("CH",2), "TH",
                     rep("CH",2), "BL",
                     rep("CH",2), "FR",
                     rep("CH",2), "AI",
                     rep("CH",2), "CM",
                     rep("CH",2), "BI",
                     rep("CH",2), "CI",
                     rep("CH",2), "DI",
                     rep("CH",2), "GF",
                     rep("CH",2), "SF",
                     rep("CH",2), "AH",
                     rep("CH",2), "EU",
                     rep("CH",2), "EI",
                     rep("CH",2), "RY",
                     rep("CH",2), "ZA",
                     rep("CH",2), "JA",
                     rep("CH",2), "HD"),
                   ncol=3, byrow=T)

group_CM <- matrix(c(rep("CM",2), "TH",
                     rep("CM",2), "BL",
                     rep("CM",2), "CH",
                     rep("CM",2), "FR",
                     rep("CM",2), "AI",
                     rep("CM",2), "BI",
                     rep("CM",2), "CI",
                     rep("CM",2), "DI",
                     rep("CM",2), "GF",
                     rep("CM",2), "SF",
                     rep("CM",2), "AH",
                     rep("CM",2), "EU",
                     rep("CM",2), "EI",
                     rep("CM",2), "RY",
                     rep("CM",2), "ZA",
                     rep("CM",2), "JA",
                     rep("CM",2), "HD"),
                   ncol=3, byrow=T)

group_BI <- matrix(c(rep("BI",2), "TH",
                     rep("BI",2), "BL",
                     rep("BI",2), "CH",
                     rep("BI",2), "FR",
                     rep("BI",2), "AI",
                     rep("BI",2), "CM",
                     rep("BI",2), "CI",
                     rep("BI",2), "DI",
                     rep("BI",2), "GF",
                     rep("BI",2), "SF",
                     rep("BI",2), "AH",
                     rep("BI",2), "EU",
                     rep("BI",2), "EI",
                     rep("BI",2), "RY",
                     rep("BI",2), "ZA",
                     rep("BI",2), "JA",
                     rep("BI",2), "HD"),
                   ncol=3, byrow=T)

group_CI <- matrix(c(rep("CI",2), "TH",
                     rep("CI",2), "BL",
                     rep("CI",2), "CH",
                     rep("CI",2), "FR",
                     rep("CI",2), "AI",
                     rep("CI",2), "CM",
                     rep("CI",2), "BI",
                     rep("CI",2), "DI",
                     rep("CI",2), "GF",
                     rep("CI",2), "SF",
                     rep("CI",2), "AH",
                     rep("CI",2), "EU",
                     rep("CI",2), "EI",
                     rep("CI",2), "RY",
                     rep("CI",2), "ZA",
                     rep("CI",2), "JA",
                     rep("CI",2), "HD"),
                   ncol=3, byrow=T)

group_DI <- matrix(c(rep("DI",2), "TH",
                     rep("DI",2), "BL",
                     rep("DI",2), "CH",
                     rep("DI",2), "FR",
                     rep("DI",2), "AI",
                     rep("DI",2), "CM",
                     rep("DI",2), "BI",
                     rep("DI",2), "CI",
                     rep("DI",2), "GF",
                     rep("DI",2), "SF",
                     rep("DI",2), "AH",
                     rep("DI",2), "EU",
                     rep("DI",2), "EI",
                     rep("DI",2), "RY",
                     rep("DI",2), "ZA",
                     rep("DI",2), "JA",
                     rep("DI",2), "HD"),
                   ncol=3, byrow=T)

group_GF <- matrix(c(rep("GF",2), "TH",
                     rep("GF",2), "BL",
                     rep("GF",2), "CH",
                     rep("GF",2), "FR",
                     rep("GF",2), "AI",
                     rep("GF",2), "CM",
                     rep("GF",2), "BI",
                     rep("GF",2), "CI",
                     rep("GF",2), "DI",
                     rep("GF",2), "SF",
                     rep("GF",2), "AH",
                     rep("GF",2), "EU",
                     rep("GF",2), "EI",
                     rep("GF",2), "RY",
                     rep("GF",2), "ZA",
                     rep("GF",2), "JA",
                     rep("GF",2), "HD"),
                   ncol=3, byrow=T)

group_SF <- matrix(c(rep("SF",2), "TH",
                     rep("SF",2), "BL",
                     rep("SF",2), "CH",
                     rep("SF",2), "FR",
                     rep("SF",2), "AI",
                     rep("SF",2), "CM",
                     rep("SF",2), "BI",
                     rep("SF",2), "CI",
                     rep("SF",2), "DI",
                     rep("SF",2), "GF",
                     rep("SF",2), "AH",
                     rep("SF",2), "EU",
                     rep("SF",2), "EI",
                     rep("SF",2), "RY",
                     rep("SF",2), "ZA",
                     rep("SF",2), "JA",
                     rep("SF",2), "HD"),
                   ncol=3, byrow=T)

group_EU <- matrix(c(rep("EU",2), "TH",
                     rep("EU",2), "BL",
                     rep("EU",2), "CH",
                     rep("EU",2), "FR",
                     rep("EU",2), "AI",
                     rep("EU",2), "CM",
                     rep("EU",2), "BI",
                     rep("EU",2), "CI",
                     rep("EU",2), "DI",
                     rep("EU",2), "GF",
                     rep("EU",2), "SF",
                     rep("EU",2), "AH",
                     rep("EU",2), "EI",
                     rep("EU",2), "RY",
                     rep("EU",2), "ZA",
                     rep("EU",2), "JA",
                     rep("EU",2), "HD"),
                   ncol=3, byrow=T)

group_EI <- matrix(c(rep("EI",2), "TH",
                     rep("EI",2), "BL",
                     rep("EI",2), "CH",
                     rep("EI",2), "FR",
                     rep("EI",2), "AI",
                     rep("EI",2), "CM",
                     rep("EI",2), "BI",
                     rep("EI",2), "CI",
                     rep("EI",2), "DI",
                     rep("EI",2), "GF",
                     rep("EI",2), "SF",
                     rep("EI",2), "AH",
                     rep("EI",2), "EU",
                     rep("EI",2), "RY",
                     rep("EI",2), "ZA",
                     rep("EI",2), "JA",
                     rep("EI",2), "HD"),
                   ncol=3, byrow=T)

group_ZA <- matrix(c(rep("ZA",2), "TH",
                     rep("ZA",2), "BL",
                     rep("ZA",2), "CH",
                     rep("ZA",2), "FR",
                     rep("ZA",2), "AI",
                     rep("ZA",2), "CM",
                     rep("ZA",2), "BI",
                     rep("ZA",2), "CI",
                     rep("ZA",2), "DI",
                     rep("ZA",2), "GF",
                     rep("ZA",2), "SF",
                     rep("ZA",2), "AH",
                     rep("ZA",2), "EU",
                     rep("ZA",2), "EI",
                     rep("ZA",2), "RY",
                     rep("ZA",2), "JA",
                     rep("ZA",2), "HD"),
                   ncol=3, byrow=T)

group_JA <- matrix(c(rep("JA",2), "TH",
                     rep("JA",2), "BL",
                     rep("JA",2), "CH",
                     rep("JA",2), "FR",
                     rep("JA",2), "AI",
                     rep("JA",2), "CM",
                     rep("JA",2), "BI",
                     rep("JA",2), "CI",
                     rep("JA",2), "DI",
                     rep("JA",2), "GF",
                     rep("JA",2), "SF",
                     rep("JA",2), "AH",
                     rep("JA",2), "EU",
                     rep("JA",2), "EI",
                     rep("JA",2), "RY",
                     rep("JA",2), "ZA",
                     rep("JA",2), "HD"),
                   ncol=3, byrow=T)

group_HD <- matrix(c(rep("HD",2), "TH",
                     rep("HD",2), "BL",
                     rep("HD",2), "CH",
                     rep("HD",2), "FR",
                     rep("HD",2), "AI",
                     rep("HD",2), "CM",
                     rep("HD",2), "BI",
                     rep("HD",2), "CI",
                     rep("HD",2), "DI",
                     rep("HD",2), "GF",
                     rep("HD",2), "SF",
                     rep("HD",2), "AH",
                     rep("HD",2), "EU",
                     rep("HD",2), "EI",
                     rep("HD",2), "RY",
                     rep("HD",2), "ZA",
                     rep("HD",2), "JA"),
                   ncol=3, byrow=T)



## Groups by Region
# group_CaMW<-matrix(c(rep("CaMW",2), "CaIN",
#                      rep("CaMW",2), "CaIS",
#                      rep("CaMW",2), "CaME",
#                      rep("CaMW",2), "GrEN",
#                      rep("CaMW",2), "GrES",
#                      rep("CaMW",2), "GrN",
#                      rep("CaMW",2), "Sib"),
#                    ncol=3, byrow=T)
# group_CaME<-matrix(c(rep("CaME",2), "CaIN",
#                      rep("CaME",2), "CaIS",
#                      rep("CaME",2), "CaMW",
#                      rep("CaME",2), "GrEN",
#                      rep("CaME",2), "GrES",
#                      rep("CaME",2), "GrN",
#                      rep("CaME",2), "Sib"),
#                    ncol=3, byrow=T)
# group_CaIS<-matrix(c(rep("CaIS",2), "CaIN",
#                      rep("CaIS",2), "CaME",
#                      rep("CaIS",2), "CaMW",
#                      rep("CaIS",2), "GrEN",
#                      rep("CaIS",2), "GrES",
#                      rep("CaIS",2), "GrN",
#                      rep("CaIS",2), "Sib"),
#                    ncol=3, byrow=T)
# group_CaIN<-matrix(c(rep("CaIN",2), "CaIS",
#                      rep("CaIN",2), "CaME",
#                      rep("CaIN",2), "CaMW",
#                      rep("CaIN",2), "GrEN",
#                      rep("CaIN",2), "GrES",
#                      rep("CaIN",2), "GrN",
#                    rep("CaIN",2), "Sib"),
#                  ncol=3, byrow=T)
# group_GrEN<-matrix(c(rep("GrEN",2), "CaIN",
#                      rep("GrEN",2), "CaIS",
#                      rep("GrEN",2), "CaME",
#                      rep("GrEN",2), "CaMW",
#                      rep("GrEN",2), "GrES",
#                      rep("GrEN",2), "GrN",
#                    rep("GrEN",2), "Sib"),
#                    ncol=3, byrow=T)
# group_GrES<-matrix(c(rep("GrES",2), "CaIN",
#                      rep("GrES",2), "CaIS",
#                      rep("GrES",2), "CaME",
#                      rep("GrES",2), "CaMW",
#                      rep("GrES",2), "GrEN",
#                      rep("GrES",2), "GrN",
#                    rep("GrES",2), "Sib"),
#                    ncol=3, byrow=T)



groups<-ls(pattern = "group_")


## CHOOSE COLOR SCHEME BASED ON GROUPING
muskoxRegCol <- c("CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853",
                   "GrN" = "#7B5D50","GrEN" = "#D02A7B", "GrES" = "#F17A30")

muskoxLocCol <- c("KU" = "#673ab7","TH" ="mediumblue","BL" ="royalblue3","CH" = "#5d9ddc","FR" = "#006064",
                  "AI" = "#009688", "CM"= "#4db6ac","BI"= "paleturquoise3","CI" ="#9bdcdd",
                  "DI"  ="darkseagreen4","GF"= "darkseagreen3","SF"= "darkolivegreen3", "AH"= "palegreen3",
                  "EU" =  "mediumseagreen", "EI" = "seagreen","RY"  ="indianred4","ZA"  ="deeppink3",
                  "JA" = "coral","HD"  ="tomato2")

## PLOT GROUP COMPARISONS: CHANGE MATRIX AND NROW TO MATCH THE NUMBER OF GROUPS!! (PNG OR PDF, YOUR CHOICE)
#png("poplike_all_multipanel3.png",width=400, height = 800)
pdf("multipanel_muscox_location.pdf",width=6, height = 20)
layout(mat=matrix(c(1,2,3,4,5,6),nrow = 6, ncol=1))
for (group in groups) {
  dstatFixH1H2ChangeH3=data.frame()
  group=get(group)
  dstatFixH1H2ChangeH3 <- do.call(rbind,apply(group, 1, getDstat3pop, dstatdf=dstatdf))
  dstatFixH1H2ChangeH3$poptree <- factor(dstatFixH1H2ChangeH3$poptree, levels=apply(group,1,paste, collapse="."))
  lalala<-chartr('.', ',',dstatFixH1H2ChangeH3$poptree)
  lalala<-sapply(lalala, gsub, pattern = ",", replacement = ", ", fixed = TRUE)
  dstatFixH1H2ChangeH3<-cbind(dstatFixH1H2ChangeH3,lalala)
  #The last population name only
  hmmm=strsplit(as.character(dstatFixH1H2ChangeH3$poptree), split = ".",fixed = TRUE)
  hehehe<-data.frame()
  hehehe[1,1]<-"yay"
  for (i in 1:length(hmmm)) {
    blarb<-hmmm[[i]][3]
    hehehe<-rbind(hehehe,blarb)
  }
  ah=length(hmmm)+1
  hehehe<-as.data.frame(hehehe[2:ah,])
  names(hehehe)[1]<-"yay"
  dstatFixH1H2ChangeH3<-cbind(dstatFixH1H2ChangeH3,hehehe)
  
  leg_scale=max(c(abs(min(dstatsSibOut$Z)),abs(max(dstatsSibOut$Z))))
  #if (leg_scale > 3 ) leg_scale=leg_scale else leg_scale=4
  
  
  par(mar=c(5,18,4,4))
  boxplot(dstatFixH1H2ChangeH3$Z ~ dstatFixH1H2ChangeH3$yay,
          boxwex=0.5,
          main=group[1],
          horizontal=T,las=2,ylim=c(-leg_scale,leg_scale),ylab="", xlab=expression(italic(Z)*"-score"),cex.axis=1.2,cex.lab=1.5,
          #names=apply(as.matrix(strsplit(as.character(dstatFixH1H2ChangeH3$poptree),split = ".",fixed = TRUE)),1, paste, collapse=", "),
          col=alpha(muskoxRegCol[group[1]], alpha=50),
          border=muskoxRegCol[group[1]])
  abline(v=c(-3,3),lty=2)
}
dev.off()
dev.off()


#png("poplike_all_multipanel3.png",width=400, height = 800)
pdf("multipanel_muscox_location.pdf",width=12, height = 20)
layout(mat=matrix(c(12,14,2,3,1,4,5,6,9,13,8,7,15,11,10),nrow = 5, ncol=3))
#layout(mat=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),nrow = 5, ncol=3))
for (group in groups) {
  dstatFixH1H2ChangeH3=data.frame()
  group=get(group)
  dstatFixH1H2ChangeH3 <- do.call(rbind,apply(group, 1, getDstat3pop, dstatdf=dstatdf))
  dstatFixH1H2ChangeH3$poptree <- factor(dstatFixH1H2ChangeH3$poptree, levels=apply(group,1,paste, collapse="."))
  lalala<-chartr('.', ',',dstatFixH1H2ChangeH3$poptree)
  lalala<-sapply(lalala, gsub, pattern = ",", replacement = ", ", fixed = TRUE)
  dstatFixH1H2ChangeH3<-cbind(dstatFixH1H2ChangeH3,lalala)
  #The last population name only
  hmmm=strsplit(as.character(dstatFixH1H2ChangeH3$poptree), split = ".",fixed = TRUE)
  hehehe<-data.frame()
  hehehe[1,1]<-"yay"
  for (i in 1:length(hmmm)) {
    blarb<-hmmm[[i]][3]
    hehehe<-rbind(hehehe,blarb)
  }
  ah=length(hmmm)+1
  hehehe<-as.data.frame(hehehe[2:ah,])
  names(hehehe)[1]<-"yay"
  dstatFixH1H2ChangeH3<-cbind(dstatFixH1H2ChangeH3,hehehe)
  
  leg_scale=leg_scale=15
    #leg_scale=max(c(abs(min(dstatsSibOut$Z)),abs(max(dstatsSibOut$Z))))
  #if (leg_scale > 3 ) leg_scale=leg_scale else leg_scale=4
  
  
  par(mar=c(3,3,4,4))
  boxplot(dstatFixH1H2ChangeH3$Z ~ dstatFixH1H2ChangeH3$yay,
          boxwex=0.5,
          main=group[1],
          horizontal=T,las=2,ylim=c(-leg_scale,leg_scale),ylab="", xlab=expression(italic(Z)*"-score"),cex.axis=1.2,cex.lab=1.5,
          #names=apply(as.matrix(strsplit(as.character(dstatFixH1H2ChangeH3$poptree),split = ".",fixed = TRUE)),1, paste, collapse=", "),
          col=alpha(muskoxLocCol[group[1]], alpha=50),
          border=muskoxLocCol[group[1]])
  abline(v=c(-3,3),lty=2)
}
dev.off()
dev.off()
