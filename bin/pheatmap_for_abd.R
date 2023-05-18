setwd("G:\\project\\YM1XLD211101MV-西沙永乐龙洞生态系统结构和生物无氧适应机制\\virus\\5.abundance\\ylbh_scs\\heatmap")
library(pheatmap)

source("utils.R")
source("color_set.R")

df<-read.table(
  "ylbh_viral_sum_abd_taxa.m.txt",header=T,sep="\t",
  row.names=1,check.names=F,quote=""
)

samp_goup_df<-read.table(
  "samp.info.xls",header=T,sep="\t",
  check.names=F,quote=""
)

colnames(samp_goup_df)<-c("sample","group")
SampNames<-samp_goup_df$sample
ColNames<-samp_goup_df$group
RowNames<-rownames(df)

ann_colors <- list(
  group=generate_named_vector(unique(ColNames),grp_color)
)
annotation_col<-data.frame(Group=factor(ColNames))
if('Source' %in% colnames(df)){
  Uniq_Source=unique(df$Source)"#000000"
  annotation_row<-data.frame(Source=factor(df$Source))
  ann_colors <- list(ann_colors,
    Source=generate_named_vector(Uniq_Source,src_color)
  )
  rownames(annotation_row)<-RowNames
}else{
  annotation_row<-NA
}

df<-subset(df,select=SampNames)
rownames(annotation_col)=names(df)

plot<-pheatmap(
  log10(df+1),
  cluster_row=TRUE,
  cluster_col=FALSE,
  show_rownames=F,
  #scale = "row",#参数归一化
  #clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward','ward.D','ward.D2','single','complete','average','mcquitty','median' or 'centroid')
  clustering_method="complete",
  #clustering_distance_rows="correlation"参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
  clustering_distance_rows="correlation",
  #color = colorRampPalette(c("blue","white","red"))(50),
  #color = colorRampPalette(c("navy","white","firebrick3"))(50),
  #color = colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100),
  annotation_colors=ann_colors,
  annotation_col=annotation_col,
  annotation_row=annotation_row
  #display_numbers=matrix(ifelse(dat > 0.01,"*",""),nrow(dat))
)
plot