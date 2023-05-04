#!/usr/bin/env Rscript
#coding:utf8
#yangm@idsse.ac.cn

if (!require("tools")) {install.packages("tools")}
if (!require("reshape2")) {install.packages("reshape2")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("RColorBrewer")) {install.packages("RColorBrewer")}

library(tools)
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
#color1 <-c(
#    "#9ACD32","#EECFA1","#4F94CD","#8B636C","#EE9A00",
#    "#FF00FF","#EECFA1","#0000FF","#BF3EFF","#00FFFF",
#    "#A52A2A","#FFFF00","#FF00FF","#00FF00","#BEBEBE",
#    "#FF0000","#EE9A00","#006400","#556B2F","#556B2F",
#    "#CD950C","#00688B","#8B795E","#458B74","#9ACDF2"
#)

#color2 <- c(
#    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
#    "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
#    "#ccebc5", "#ffed6f", "#8C9EFF", "#9b9b9b", "#ffffe5",
#    "#fdb462", "#bebada", "#fb8072", "#80b1d3", "#ffffb3",
#    "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
#--------------
#    "#9ACD32","#EECFA1","#4F94CD","#8B636C","#EE9A00",
#    "#FF00FF","#EECFA1","#0000FF","#BF3EFF","#00FFFF",
#    "#A52A2A","#FFFF00","#FF00FF","#00FF00","#BEBEBE",
#    "#FF0000","#EE9A00","#006400","#556B2F","#556B2F",
#    "#CD950C","#00688B","#8B795E","#458B74","#9ACDF2"
#)
color <- c("#8dd3c7", "#fdb462", "#b3de69", "#82c1a8", "#d9d9d9",
           "#9ac73f", "#c34b48", "#a163ad", "#ff1493", "#00c1ff",
           "#a87c00", "#619cff", "#f8766d", "#008000", "#964500",
           "#f6b3c3", "#5b8cb2", "#d2d28e", "#b297c9", "#ccebc5",
           "#cc3a21", "#00ba38", "#a15d98", "#004a7f", "#ffc0cb",
           "#848482", "#00ffff", "#ffed6f", "#9b9b9b", "#ffffe5",
           "#bebada", "#fb8072", "#e6a141", "#80b1d3", "#f63200",
           "#fccde5", "#bc80bd", "#5da5b3", "#ff00ff", "#556b2f", 
           "#008080", "#b41e00", "#6c7c32", "#cc8899", "#8f00ff",
           "#ffa500", "#3b3b3b", "#ff0000", "#006400", "#00a3a3",
           "#cd950c", "#00688b", "#8b795e", "#458b74", "#9acdf2",
           "#ee9a00", "#00ff00", "#bebebe", "#fffac8", "#beaed4",
           "#cd950c", "#00688B")

plot<-ggplot(df_p,aes(fill=Type,y=Relative_Abundance,x=SampleID))+
  geom_bar(stat="identity",position="fill",width=.8,color="black")+
  scale_fill_manual(values=color)+
  theme_bw()+
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),
        axis.title=element_text(size=10),
        axis.text=element_text(size=10))
outprefix=file_path_sans_ext(basename(argv[1]))
pdf(paste(argv[2],"/",outprefix,".barplot.pdf",sep=""),width=10,height=8)
plot
dev.off()
