######## PCA
require(ggplot2)
require(tidyverse)
require(ggrepel)
require(data.table)

## Upload covariance matrix
muskox_covariance <- as.matrix(read.table("pcangsd_wSib.cov"))

## Upload pop info
muskoxSum1 <- read.csv("muskox_bamlist.txt", sep = "\t", header = TRUE)
muskoxSum2 <- read.csv("samples_wSib.txt", sep = "\t", header = TRUE)
muskoxSum <- inner_join(muskoxSum2,muskoxSum1,by="Sample")

## Apply eigen function
e <- eigen(muskox_covariance)

## Calculate the proportion of variance in the PCs
varPC1 = format(round(e$values[1]/sum(e$values)*100, 2), nsmall = 2)
varPC2 = format(round(e$values[2]/sum(e$values)*100, 2), nsmall = 2)
varPC3 = format(round(e$values[3]/sum(e$values)*100, 2), nsmall = 2)
varPC4 = format(round(e$values[4]/sum(e$values)*100, 2), nsmall = 2)

## Convert vectors as data frame
edf <- as.data.frame(e$vectors)

## Combine with info on the samples
edf <- cbind(muskoxSum$Sample, muskoxSum$Location, muskoxSum$Region, edf)
names(edf)[names(edf) == 'muskoxSum$Sample'] <- 'Sample'
names(edf)[names(edf) == 'muskoxSum$Location'] <- 'Location'
names(edf)[names(edf) == 'muskoxSum$Region'] <- 'Region'

## Define location levels
edf$Location <- factor(edf$Location, levels = c("SIB","KU","TH","BL","CH","FR","AI","CM","BI","CI","DI","GF","SF","AH",
                                                "EU","EI","RY","ZA","JA","HD"))
edf$Region <- factor(edf$Region, levels = c("SIB","CaMW", "CaME","CaIS","CaIN","GrN", "GrEN","GrES"))

## Define color palette
muskoxLocCol <- c("#944ead","#673ab7","mediumblue","royalblue3", "#5d9ddc","#006064",
                  "#009688","#4db6ac","paleturquoise3","#9bdcdd",
                  "darkseagreen4","darkseagreen3","darkolivegreen3","palegreen3",
                  "mediumseagreen","seagreen","indianred4","deeppink3",
                  "coral","tomato2")
muskoxRegCol <- c("SIB" = "#967bb6","CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853",
                  "GrN" = "#7B5D50", "GrEN" = "#D02A7B", "GrES" = "#F17A30")

## Plot PC1 vs PC2
PCA1v2_pop <- ggplot(edf, aes(x=V1, y=V2)) + geom_point(aes(x=V1, y=V2, fill = edf$Location),shape = 21, size =10) + 
        scale_fill_manual(name = "Populations", labels = levels(edf$Location), breaks = levels(edf$Location), values = muskoxLocCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC1 (",varPC1,"%)",sep=""), y=paste("PC2 (",varPC2,"%)",sep="")) + 
#       geom_text_repel(label = edf$Sample, size = 7) + 
       theme(axis.text=element_text(size = rel(2.2)), axis.title = element_text(size =rel(2.5)), axis.title.x = element_text(vjust = -1),
              legend.text = element_text(size = rel(2.2)), legend.title = element_text(size = rel(2.2)))
PCA1v2_pop

PCA2v3_pop <- ggplot(edf, aes(x=V2, y=V3)) + geom_point(aes(x=V2, y=V3, fill = edf$Location),shape = 21, size =10) + 
        scale_fill_manual(name = "Populations", labels = levels(edf$Location), breaks = levels(edf$Location), values = muskoxLocCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC2 (",varPC2,"%)",sep=""), y=paste("PC3 (",varPC3,"%)",sep="")) + 
  #     geom_text_repel(label = edf$Sample, size = 7) + 
        theme(axis.text=element_text(size = rel(2.2)), axis.title = element_text(size =rel(2.5)), axis.title.x = element_text(vjust = -1),
              legend.text = element_text(size = rel(2.2)), legend.title = element_text(size = rel(2.2)))
PCA2v3_pop

PCA3v4_pop <- ggplot(edf, aes(x=V3, y=V4)) + geom_point(aes(x=V3, y=V4, fill = edf$Location),shape = 21, size =10) + 
        scale_fill_manual(name = "Populations", labels = levels(edf$Location), breaks = levels(edf$Location), values = muskoxLocCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC3 (",varPC3,"%)",sep=""), y=paste("PC4 (",varPC4,"%)",sep="")) + 
#        geom_text_repel(label = edf$Sample, size = 7) + 
        theme(axis.text=element_text(size = rel(2.2)), axis.title = element_text(size =rel(2.5)), axis.title.x = element_text(vjust = -1),
              legend.text = element_text(size = rel(2.2)), legend.title = element_text(size = rel(2.2)))
PCA3v4_pop

## Plot PCs
PCA1v2_reg <- ggplot(edf, aes(x=V1, y=V2)) + geom_point(aes(x=V1, y=V2, fill = edf$Region),shape = 21, size =10) + 
        scale_color_manual(name = "Populations", labels = levels(edf$Region), breaks = levels(edf$Region), values = muskoxRegCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC1 (",varPC1,"%)",sep=""), y=paste("PC2 (",varPC2,"%)",sep="")) + 
        theme(axis.text=element_text(size = rel(2)), axis.title = element_text(size =rel(2)), axis.title.x = element_text(vjust = -1),
          legend.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2))) 
PCA1v2_reg

PCA2v3_reg <- ggplot(edf, aes(x=V2, y=V3)) + geom_point(aes(x=V2, y=V3, fill = edf$Region),shape = 21, size =10) + 
        scale_color_manual(name = "Populations", labels = levels(edf$Region), breaks = levels(edf$Region), values = muskoxRegCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC2 (",varPC2,"%)",sep=""), y=paste("PC3 (",varPC3,"%)",sep="")) + 
        theme(axis.text=element_text(size = rel(2)), axis.title = element_text(size =rel(2)), axis.title.x = element_text(vjust = -1),
              legend.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2))) 
PCA2v3_reg

PCA3v4_reg <- ggplot(edf, aes(x=V3, y=V4)) + geom_point(aes(x=V3, y=V4, fill = edf$Region),shape = 21, size =10) + 
        scale_color_manual(name = "Populations", labels = levels(edf$Region), breaks = levels(edf$Region), values = muskoxRegCol, aesthetics = c("colour","fill")) + 
        theme_classic() + labs(x=paste("PC3 (",varPC3,"%)",sep=""), y=paste("PC4 (",varPC4,"%)",sep="")) + 
        theme(axis.text=element_text(size = rel(2)), axis.title = element_text(size =rel(2)), axis.title.x = element_text(vjust = -1),
              legend.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2))) 
PCA3v4_reg
