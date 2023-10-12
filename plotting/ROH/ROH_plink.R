## load modules
require(ggplot2)
require(ggrepel)
require(plyr)
require(dplyr)
require(forcats)
require(tidyverse)
require(scales)

## Define color palette
muskoxLocCol <- c("#967bb6","#673ab7","mediumblue","royalblue3", "#5d9ddc",
                  "#009688","#4db6ac","paleturquoise3","#9bdcdd","#006064",
                  "darkseagreen4","darkseagreen3","darkolivegreen3","palegreen3",
                  "mediumseagreen","seagreen","indianred4","deeppink3",
                  "coral","tomato2")
muskoxRegCol <- c("SIB"="#967bb6","CaMW" = "#8444DE", "CaME" = "#618DFF", "CaIS" = "#41A99C", "CaIN" = "#368853",
                  "GrN" = "#7B5D50","GrEN" = "#D02A7B", "GrES" = "#F17A30")

#### PLINK

roh <- read.table("roh.txt", 
                   header=TRUE)

######################### DEFINE POPULATIONS

muskoxLocCol <- c("#967bb6","#673ab7","mediumblue","royalblue3", "#5d9ddc",
                  "#009688","#4db6ac","paleturquoise3","#9bdcdd",
                  "darkolivegreen3","palegreen3",
                  "mediumseagreen","seagreen","indianred4","deeppink3",
                  "coral","tomato2")
## order levels
roh$Location <- factor(roh$Location, levels = c("SIB","KU","TH","BL","CH","AI","CM","BI","CI","SF","AH",
                                                "EU","EI","RY","ZA","JA","HD"))
roh$Region <- factor(roh$Region, levels = c("SIB","CaMW", "CaME","CaIS","CaIN","GrN","GrEN","GrES"))

########################### PLOTTING

## plot

rohPlot <- ggplot(roh, aes(x = as.factor(roh$Location), y= as.numeric(roh$PLINK_Froh), fill = roh$Location)) + 
  geom_bar(aes(x = as.factor(roh$Location)), size = 4, stat="summary", fun = "mean") + ylab(expression(F[ROH])) + xlab("") +
  stat_summary(fun.data = "mean_se", geom = "errorbar",width = .2) + 
  theme(axis.text.x = element_text(size = rel(2.6), colour = "grey20", angle = 45, vjust = 0.5, hjust = 0.5), plot.margin=unit(c(1.2,1,1,1),"cm"),
        axis.text.y = element_text(colour = "black", size = rel(2.5)), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.y  = element_text(size = rel(2.7), margin = margin(t = 0, r = 40, b = 0, l = 0)), axis.title.x = element_blank(),
        #        plot.title = element_text(size = rel(3), face = "bold"), 
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_blank(),
        #        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.title = element_text(size = rel(2.5)), legend.text = element_text(size = rel(2.2))) +
  scale_y_continuous(limits = c(0.00,1), expand = c(0,0), breaks = seq(0,8,0.25)) +
  #  geom_text_repel(label = het$Sample, size = 6, point.padding = 1)
  scale_fill_manual(name = "Locations", values = muskoxLocCol, labels = levels(roh$Location))
rohPlot


### plot stacked
## plot locations
roh <- read.table("plink/PLINKsummary.txt", 
                  sep = "", header=TRUE)
roh$Location <- factor(roh$Location, levels = c("SIB","KU","TH","BL","CH","AI","CM","BI","CI","FR","DI","GF","SF","AH",
                                                "EU","EI","RY","ZA","JA","HD"))
roh$Region <- factor(roh$Region, levels = c("SIB","CaMW", "CaME","CaIS","CaIN","GrN","GrEN","GrES"))
roh$Size <- factor(roh$Size, levels = c("0.5-1mb","1-2mb","2-5mb","5-10mb","≥10mb"))
roh$Sample <- factor(roh$Sample, levels = c("EL012", "KU-11_64","KU-12_70","KU-13_38","KU-21_43",
                                            "KU-23_3","KU-27_2","KU-2_57","KU-31_46","KU-32_48",
                                            "KU-39_56","KU-4_72","TH-27_69","BL-10_175","BL-11_9","BL-13_8",
                                            "BL-15_31","BL-17_11","BL-1_7","BL-20_10","BL-23_32",
                                            "BL-24_33","BL-2_26","BL-4_28","BL-6_29","CH-3_62",
                                            "CH-7_45","AI-1_119","CM-1_61","CM-2_63","BI-1_162",
                                            "BI-26_18","BI-28_99","BI-32_19","BI-35_20","BI-36_21",
                                            "CI-8_14","SF-1_42","SF-3_49","AH-1_40","EU-7_67",
                                            "EI-23_12","EI-28_17", "RY-3","MO25","MO28","MO29",
                                            "MO32","MO33","MO34","MO35","MO42","MO46",
                                            "MO49", "MO50","JA-1_34","JA-4_36","HD-1","HD-2","HD-3",
                                            "HD-5"))

cbf_1 <- c("#009E73", "#56B4E9", "#999999", 
           "#E69F00", "#CC79A7")

rohPlot <- ggplot(roh, aes(x = as.factor(roh$Sample), y= as.numeric(roh$ROH/1000), fill = as.factor(roh$Size))) + 
  geom_bar(position = "stack", stat="identity") + theme(axis.text.x = element_text(angle = 90, size = rel(1.7), vjust = 0.5), 
                                                        axis.title.y  = element_text(size = rel(2)),
                                                        axis.text.y  = element_text(size = rel(2)),
                                                        axis.title.x = element_text(size = rel(2)),
                                                        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.7))) +
                                    xlab("Samples") + ylab("Size of genome in ROH in MB") +
  scale_fill_manual(name = "ROH size category",values = cbf_1, labels = levels(roh$Size))
rohPlot
