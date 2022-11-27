#!/usr/bin/env Rscript
#coding:utf8

if (!require("reshape2")) {install.packages("reshape2")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("RColorBrewer")) {install.packages("RColorBrewer")}
library(reshape2)
library(ggplot2)
library(RColorBrewer)

argv<-commandArgs(T)

df<-read.table(
    argv[1],header=TRUE,sep="\t",
    row.names=1,check.names=F,quote=""
)

df_t<-as.data.frame(t(df))
df_t<-cbind("SampleID"=(rownames(df_t)),df_t)
rownames(df_t)<-NULL
df_p<-melt(
    df_t,id.vars="SampleID",variable.name="Type",
    value.name="Relative_Abundance"
)
df_p$SampleID<-factor(
    df_p$SampleID,
    levels=df_t$SampleID
)
#"#556B2F","#556B2F","#CD950C","#00688B","#8B795E","#458B74",
color<-c(
           "#9ACD32","#EECFA1","#4F94CD","#8B636C","#EE9A00",
           "#FF00FF","#EECFA1","#0000FF","#BF3EFF","#00FFFF",
           "#A52A2A","#FFFF00","#FF00FF","#00FF00","#BEBEBE",
           "#FF0000","#EE9A00","#006400"
)

plot<-ggplot(df_p,aes(fill=Type,y=Relative_Abundance,x=SampleID))+
  geom_bar(stat="identity",position="fill",width =.8,color="black") +
  scale_fill_manual(values =color)+
  theme_bw()+
  theme(legend.text =element_text(size=10),
        legend.title =element_text(size=10),
        axis.text.x =element_text(angle=90,hjust=0.5,vjust=0.5),
        axis.title =element_text(size=10),
        axis.text =element_text(size=10))
pdf(paste(argv[2],"/taxa_relative_abundance_baplot.pdf",sep=""),width=10,height=8)
plot
