#!/usr/bin/env Rscript
#coding:utf8

if (!require("reshape2")) {install.packages("reshape2")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("RColorBrewer")) {install.packages("RColorBrewer")}
library(reshape2)
library(ggplot2)
library(RColorBrewer)

argv <- commandArgs(T)

dat <- read.table(
    argv[1], header=TRUE, sep="\t", 
    row.names=1, check.names=F, quote=""
)

dat_a <- dat[seq(1:18),]
dat_b <- as.data.frame(t(dat_a))
dat_b <- cbind("SampleID"=(rownames(dat_b)),dat_b)
rownames(dat_b) <- NULL
dat_c <- melt(dat_b, id.vars="SampleID", variable.name="Type", value.name="Relative_Abundance")
dat_c$SampleID <- factor(
    dat_c$SampleID,
    levels=c("V1","V2","V3",
           "CB148MP1_CB","CB148MP2_CB","CB200MS1_CB","CB237MS1_CB",
           "CB247MS1_CB","CB267MS1_CB","CB267MS2_CB","CB900MS1_CB",
           "BH0560M_AH","BH0960M_AH","BH0995M_AH","BH05106M_AH"
    )
)
#"#556B2F","#556B2F","#CD950C","#00688B","#8B795E","#458B74",
color <- c(
           "#9ACD32","#EECFA1","#4F94CD","#8B636C","#EE9A00",
           "#FF00FF","#EECFA1","#0000FF","#BF3EFF","#00FFFF",
           "#A52A2A","#FFFF00","#FF00FF","#00FF00","#BEBEBE",
           "#FF0000","#EE9A00","#006400"
)

plot <- ggplot(dat_c, aes(fill= Type, y= Relative_Abundance, x= SampleID))+
  geom_bar(stat="identity", position="fill",width = .8,color="black") +
  scale_fill_manual(values = color)+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
pdf(paste(argv[2],"/taxa_relative_abundance_baplot.pdf", width=10, height=8)
plot
