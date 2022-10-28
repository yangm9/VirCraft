#!/usr/bin/env Rscript
#coding:utf8

if (!require("pheatmap")) {install.packages("pheatmap")}
library(pheatmap)

argv <- commandArgs(T)

if (length(argv) < 3){
    stop("inputs: <merged_tpm_xls> <samp_group_xls> <heatmap_dir>")
}

df <- read.table(
    argv[1], header=TRUE, sep="\t", 
    row.names=1, check.names=F, quote=""
)

samp_goup_df <- read.table(
    argv[2], header=TRUE, sep="\t", 
    check.names=F, quote=""
)

ColNames <- samp_goup_df$GroupName
RowNames <- rownames(df)

annotation_col<-data.frame(
  Group=factor(ColNames)
)

annotation_row = data.frame(Source=factor(df$Source))
rownames(annotation_row) <- RowNames
df <- subset(df, select=-c(Order, Family, Source, Length))
rownames(annotation_col)=names(df)

plot <- pheatmap(
  log10(df+1),
  cluster_row=TRUE, 
  cluster_col=FALSE,
  show_rownames = F,
  #scale = "row", #参数归一化
  #clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid')
  clustering_method="complete",
  #clustering_distance_rows="correlation"参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
  clustering_distance_rows="correlation",
  #color = colorRampPalette(c("blue","white","red"))(50),
  #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  #color = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
  annotation_col=annotation_col,
  annotation_row=annotation_row
  #display_numbers=matrix(ifelse(dat > 0.01, "*", ""), nrow(dat))
)

pdf(paste(argv[3], '/tpm_heatmap.pdf', sep=''), width=10,height=8)
plot
dev.off()
#png(paste(argv[3], '/tpm_heatmap.png', sep=''), width=1000,height=800)
#plot
#dev.off()
