#!/usr/bin/env Rscript
#coding:utf8
#yangm@idsse.ac.cn

if (!require("tools")) {install.packages("tools")}
if (!require("vegan")) {install.packages("vegan")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("patchwork")) {install.packages("patchwork")}

library(tools)
library(vegan)
library(ggplot2)
library(patchwork)

argv<-commandArgs(T)

if (length(argv)<2){
    stop("inputs: <all_abundance.xls> <alpha_diversity.xls>")
}

otu<-read.table(argv[1],sep="\t",header=T,
                check.names=FALSE,row.names=1)
otu_t<-t(otu)
Shannon=diversity(otu_t,"shannon") 
Simpson=diversity(otu_t,"simpson") 
InSimpson=diversity(otu_t,"inv")
specN=specnumber(otu_t)
Pielou=Shannon/log(specN)  #计算Pielou均匀度指数
alpha_df=t(rbind(specN,Pielou,Shannon,Simpson,InSimpson))
alpha_df<-as.data.frame(alpha_df)
alpha_df$Sample<-rownames(alpha_df)
rownames(alpha_df)<-NULL
alpha_df<-alpha_df[,c(ncol(alpha_df),1:(ncol(alpha_df)-1))] #最后Index列移动到第一列
outprefix=file_path_sans_ext(basename(argv[1]))
alpha_tab=paste(argv[3],'/',outprefix,'.alpha_diversity.xls',sep='')
write.table(alpha_df,alpha_tab,sep='\t',row.names=TRUE)

samp_group_df<-read.table(argv[2],header=T,sep="\t", #samp.info.xls
                          check.names=F,quote="")
samp_group_df<-samp_group_df[,1:2]
alpha_df<-merge(samp_group_df,alpha_df,by.x="Sample",by.y="Sample")

color1<-c("#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
color2<-c("#e377c2", "#7f7f7f", "#bcbd22", "#9edae5")
color3<-c("#aec7e8", "#ffbb78", "#98df8a", "#ff9896")
color4<-c("#f7b6d2", "#17becf", "#dbdb8d", "#c7c7c7")

p1<-ggplot(alpha_df,aes(x=Group,y=specN)) +
    geom_boxplot(outlier.size=1,fill=color1) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill='transparent',color='black')) +
    labs(x='',y='specN')

p2<-ggplot(alpha_df,aes(x=Group,y=Shannon)) +
    geom_boxplot(outlier.size=1,fill = color1) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill='transparent',color='black')) +
    labs(x='', y ='Shannon')

p3<-ggplot(alpha_df, aes(x=Group, y=Simpson)) +
    geom_boxplot(outlier.size = 1, fill = color1) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill='transparent',color='black')) +
    labs(x='', y ='Simpson')

p4<-ggplot(alpha_df,aes(x=Group,y=Pielou)) +
    geom_boxplot(outlier.size=1,fill=color1) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill='transparent',color='black')) +
    labs(x='',y='Pielou')

pdf(paste(argv[3],'/',outprefix,'.alpha_diversity.boxplot.pdf',sep=''),width=10,height=8)
p1+p2+p3+p4
dev.off()
