## load modules
require(ggplot2)
require(ggrepel)
require(dplyr)
require(scales)

###################### THE GC CALL option

## load file with heterozygosity
het_filt <- read.table("het_comparison_HQ_2_table_comparison.txt", 
                  sep = "", header=TRUE)

het_filt$Location <- factor(het_filt$Location, levels = c("SIB","KU","TH","BL","CH","AI","CM","BI","CI","SF","AH",
                                                          "EU","EI","RY","ZA","JA","HD"))
muskoxLocCol <- c("#967bb6","#673ab7","mediumblue","royalblue3", "#5d9ddc",
                  "#009688","#4db6ac","paleturquoise3","#9bdcdd",
                  "darkolivegreen3","palegreen3",
                  "mediumseagreen","seagreen","indianred4","deeppink3",
                  "coral","tomato2")


## plot locations
hetBoxPlot <- ggplot(het_filt, aes(x = Location, y=as.numeric(het_filt$Het_muskoxmap_GC), color = het_filt$Location)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") +
  #  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), plot.margin=unit(c(1.2,1,1,1),"cm"),
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000005,0.00099), expand = c(0,0), breaks = seq(0,8,0.0001)) +
  scale_color_manual(name = "Locations", values = muskoxLocCol, breaks = levels(as.factor(het_filt$Location)), labels = levels(het_filt$Location))
hetBoxPlot

## plot locations log
hetBoxPlot <- ggplot(het_filt, aes(x = Location, y=as.numeric(het_filt$Het_muskoxmap_GC), color = het_filt$Location)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") +
  #  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), plot.margin=unit(c(1.2,1,1,1),"cm"),
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.y = element_line(colour = "grey70"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000005,0.00099), expand = c(0,0), trans='log10', breaks = c(0.000005, 0.00001, 0.00002, 0.00005,0.0001, 0.0002, 0.0004, 0.0006, 0.0008)) +
  scale_color_manual(name = "Locations", values = muskoxLocCol, breaks = levels(as.factor(het_filt$Location)), labels = levels(het_filt$Location))
hetBoxPlot

## plot locations vertical

het_filt$Location <- factor(het_filt$Location, levels = c("HD","JA","ZA","RY","EI","EU","AH","SF","CI","BI","CM","AI","CH","BL","TH","KU","SIB"))
muskoxLocCol <- c("tomato2", "coral","deeppink3","indianred4","seagreen","mediumseagreen","palegreen3","darkolivegreen3","#9bdcdd","paleturquoise3",
                  "#4db6ac","#009688","#006064","royalblue3","mediumblue","#673ab7","#967bb6")

hetBoxPlot <- ggplot(het_filt, aes(x = Location, y=as.numeric(het_filt$Het_muskoxmap_GC), color = het_filt$Location)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") + coord_flip() +
  #  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), plot.margin=unit(c(1.2,1,1,1),"cm"),
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.x = element_line(colour = "grey80"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000,0.00088), expand = c(0,0), breaks = seq(0,8,0.0001)) +
  scale_color_manual(name = "Locations", values = muskoxLocCol, breaks = levels(as.factor(het_filt$Location)), labels = levels(het_filt$Location))
hetBoxPlot

##### plot regions
het_filt$Region <- factor(het_filt$Region, levels = c("SIB","CaMW", "CaME","CaIS","CaIN","GrN","GrEN","GrES"))
muskoxRegCol <- c("SIB" = "#967bb6","CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853",
                  "GrN" = "#7B5D50","GrEN" = "#D02A7B", "GrES" = "#F17A30")

## plot regions
hetBoxPlot <- ggplot(het_filt, aes(x = het_filt$Region, y=as.numeric(het_filt$Het_muskoxmap_GC), color = het_filt$Region)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") +
  #  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), plot.margin=unit(c(1.2,1,1,1),"cm"),
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000,0.00088), expand = c(0,0), breaks = seq(0,8,0.0001)) +
  scale_color_manual(name = "Locations", values = muskoxRegCol, breaks = levels(as.factor(het_filt$Region)), labels = levels(het_filt$Region))
hetBoxPlot

########################## species comparison

het_filt <- read.table("Table_het_all_updated_morespecies_diet2.txt", 
                       sep = "\t", header=TRUE)

hetcol <- c(rep("black", 49), "#DC19A0", "#DC19A0")


ggplot(het_filt) + 
  geom_bar(mapping = aes(y = as.numeric(het_filt$Genome.wide.heterozygosity), x = reorder(as.factor(het_filt$Common.Name), as.numeric(-het_filt$Genome.wide.heterozygosity))), 
           stat = "identity", fill = c("grey70", "#EB669B")[as.factor(het_filt$Color[order(het_filt$Genome.wide.heterozygosity)])]) + 
  coord_flip() + theme_minimal() + ylab("Genome-wide heterozygosity") + xlab("") +
  theme(axis.text.x = element_text(size = rel(2.5)), axis.text.y = element_text(size = rel(2.5)), 
        axis.title.x = element_text(size = rel(2.5), margin = margin(t = 10, r = 0, b = 0, l = 0)), plot.margin=unit(c(1.2,1,1,1),"cm"))+
  scale_fill_manual(values = c("grey70", "#EB669B")[het_filt$Color[order(-het_filt$Genome.wide.heterozygosity)]])
 
  

###################### GL option

het <- read.table("all_est.ml", 
                   sep = "", header=FALSE)
het_Est <- as.data.frame(het[,3]/rowSums(het[,2:4]))
het_Names <- as.data.frame(het[,1])
het_All <- cbind(het_Names, het_Est)
colnames(het_All) <- c("Sample","Heterozygosity")

## Upload pop info and error rates
#samples <- read.csv("../muskox/remap/heterozygosity/samples.txt", sep = "\t", header = TRUE)
errorRates <- read.csv("muskox_bamlist_Sib.txt", 
                        sep = "\t", header = TRUE)

het <- inner_join(errorRates,het_All,by="Sample")

## DEFINE POPULATIONS

## order levels
het$Location <- factor(het$Location, levels = c("SIB","KU","TH","BL","CH","FR","SI","AI","CM","BI","CI","DI","GF","SF","AH",
                                                              "EU","EI","RY","BS","KC","NG","ZA","HH","YM","GI","JA","HD"))
het$Region <- factor(het$Region, levels = c("SIB","CaMW", "CaME","CaIS","CaIN","GrN","GrEN","GrES"))

## PLOTTING

muskoxLocCol <- c("#967bb6","#673ab7","mediumblue","royalblue3", "#5d9ddc","#006064",
                           "#009688","#4db6ac","paleturquoise3","#9bdcdd",
                           "darkseagreen4","darkseagreen3","darkolivegreen3","palegreen3",
                           "mediumseagreen","seagreen","indianred4","deeppink3",
                           "coral","tomato2","goldenrod2")

## plot individuals 
hetPlot <- ggplot(het, aes(x = as.factor(Sample), y= as.numeric(Heterozygosity), color = het$Location)) + 
  geom_point(aes(x = as.factor(het$Sample)), size = 4) + ylab("Heterozygosity") +
labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(0.7), colour = "grey20", angle = 90, vjust = 0.2, hjust = 0.5), 
        axis.text.y = element_text(face="bold", colour = "black", size = rel(1.8)), 
        axis.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.title = element_text(size = rel(2), face = "bold"), panel.grid.major.y = element_line(colour = "grey50"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.6))) +
  scale_y_continuous(limits = c(0.000,0.0006), expand = c(0,0), breaks = seq(0,8,0.00005)) +
#  geom_text_repel(label = het$Sample, size = 6, point.padding = 1)
  scale_color_manual(name = "Locations", values = muskoxLocCol, labels = levels(het$Location))
hetPlot

## plot individuals with Sib
hetPlotLocs <- ggplot(het, aes(x = as.factor(Sample), y=as.numeric(Heterozygosity), color = het$Location)) + 
  geom_point(aes(x = as.factor(Sample)), size = 4) +
  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(1), colour = "grey20", angle = 45, vjust = 0.2, hjust = 0.5), 
        axis.text.y = element_text(face="bold", colour = "black", size = rel(1.8)), 
        axis.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.title = element_text(size = rel(2), face = "bold"), panel.grid.major.y = element_line(colour = "grey50"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.6))) +
  scale_y_continuous(limits = c(0.000,0.005), expand = c(0,0), breaks = seq(0,8,0.0005)) +
  scale_color_manual(name = "Locations", values = muskoxLocCol, labels = levels(het$Location)) 
hetPlotLocs

## Boxplot

## plot locations
hetBoxPlot <- ggplot(het, aes(x = Location, y=as.numeric(het$Heterozygosity), color = het$Location)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") +
  #  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), 
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000,0.005), expand = c(0,0), breaks = seq(0,8,0.0005)) +
  scale_color_manual(name = "Locations", values = muskoxLocCol, breaks = levels(as.factor(het$Location)), labels = levels(het$Location))
hetBoxPlot

## plot regions
hetBoxPlot <- ggplot(het, aes(x = Region, y=as.numeric(Heterozygosity), color = Region)) + 
  geom_boxplot(lwd = 2) + theme_classic() + ylab("Heterozygosity") +
#  labs(title = "Heterozygosity", vjust = 1) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), 
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_blank(),
#        plot.title = element_text(size = rel(3), face = "bold"), 
      plot.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80"),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.000,0.005), expand = c(0,0), breaks = seq(0,8,0.0005)) +
  scale_color_manual(name = "Regions", values = muskoxRegCol, breaks = levels(as.factor(het_filt$Region)), labels = levels(het_filt$Region))
hetBoxPlot

##### CPG, transition comparison

hetcomp <- read.table("het_comparison_notrans_plot.txt", header = T, sep = "\t")
hetcomp$Location <- factor(hetcomp$Location, levels = c("SIB","KU","TH","BL","CH","FR","SI","AI","CM","BI","CI","DI","GF","SF","AH",
                                                        "EU","EI","RY","BS","KC","NG","ZA","HH","YM","GI","JA","HD"))
ggplot(hetcomp, aes(x = hetcomp$Location, y = as.numeric(hetcomp$Het), fill = hetcomp$estimate)) + geom_boxplot() +
  xlab("Locations") + ylab("Genome-wide heterozygosity (log10 scale)") + scale_fill_manual(name = "", values = c("#C33F6F","#767676","#EB669B","grey70")) + 
  scale_y_continuous(trans = 'log10')

