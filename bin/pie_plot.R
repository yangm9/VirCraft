#!/usr/bin/env Rscript
#yangm@idsse.ac.cn

argv<-commandArgs(T)

if(length(argv)<3){
    stop("inputs: <variables.xls> <Variable> <plot_dir>\n")
}

if(!require("ggplot2")){install.packages("ggplot2")}
if(!require("RColorBrewer")){install.packages("RColorBrewer")}

library(ggplot2)
library(RColorBrewer)

dt <- read.table(argv[1],header=T,sep="\t",
                 row.names=1,check.names=F,quote="")

df <- data.frame(table(dt[,argv[2]])) #Field
colnames(df) <- c("Variable", "Count")
df = df[order(df[,"Count"], decreasing = TRUE),] 
myLabel = as.vector(df[,"Variable"]) 
myLabel = paste(myLabel,
                "(", round(df[,"Count"] / sum(df[,"Count"]) * 100, 2), 
                "%)", sep = "") 
color <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
  "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
  "#ccebc5", "#ffed6f", "#8C9EFF", "#9b9b9b", "#ffffe5",
  "#fdb462", "#bebada", "#fb8072", "#80b1d3", "#ffffb3",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"
)

pie_plot <- ggplot(df,aes(x = "", y = Count, fill = Variable)) + 
  geom_bar(stat="identity", width = 1) + 
  coord_polar(theta="y") + 
  labs(x="", y="", title=paste(argv[2],"Pie Plot",sep=' ')) + 
  theme(axis.ticks=element_blank()) + 
#  theme(legend.title=element_blank(), legend.position = "top") + 
  scale_fill_manual(values=color) + 
#  geom_text(aes(y = Count/2 + c(0, cumsum(Count)[-length(Count)]), 
#                x = sum(Count)/20, label = myLabel), 
#                size = 3) + ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置 p
  theme_void() +
  theme(axis.text.x = element_blank())
pdf(paste(argv[3],'/',argv[2],'_pie_plot.pdf',sep=''),width=10,height=8)
pie_plot
dev.off()
