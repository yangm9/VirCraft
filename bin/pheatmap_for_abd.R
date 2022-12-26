#!/usr/bin/env Rscript
#为每一行添加来源或者Taxa
argv<-commandArgs(T)

if(length(argv)<3){
    stop("inputs: <merged_tpm_xls> <samp_group_xls> <heatmap_dir>\n")
}

if(!require("pheatmap")){install.packages("pheatmap")}
library(pheatmap)

df<-read.table(
    argv[1],header=T,sep="\t",
    row.names=1,check.names=F,quote=""
)

samp_goup_df<-read.table(
    argv[2],header=T,sep="\t",
    check.names=F,quote=""
)

colnames(samp_goup_df)<-c("sample","group")
SampNames<-samp_goup_df$sample
ColNames<-samp_goup_df$group
RowNames<-rownames(df)

annotation_col<-data.frame(Group=factor(ColNames))
if('Source' %in% colnames(df)){
    annotation_row<-data.frame(Source=factor(df$Source))
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
    annotation_col=annotation_col,
    annotation_row=annotation_row
    #display_numbers=matrix(ifelse(dat > 0.01,"*",""),nrow(dat))
)

pdf(paste(argv[3],'/abundance_heatmap.pdf',sep=''),width=10,height=8)
plot
dev.off()
#Sys.setenv("DISPLAY"=":0")
#png(paste(argv[3],'/tpm_heatmap.png',sep=''),width=1000,height=800)
#plot
#dev.off()
