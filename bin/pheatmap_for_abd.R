#!/usr/bin/env Rscript
#yangm@idsse.ac.cn
#为每一行添加来源或者Taxa
argv <- commandArgs(T)

if(length(argv) < 3){
    stop("inputs: <merged_tpm_tsv> <samp_group_tsv> <heatmap_dir>\n")
}

if(!require("pheatmap")){
    install.packages("pheatmap")
}

library(pheatmap)

# 生成命名向量函数
generate_named_vector <- function(keys, values) {
    named_vector <- setNames(values[1:length(keys)], keys)
    return(named_vector)
}

# 读取数据
df <- read.table(
    argv[1], header=T, sep="\t",
    row.names=1, check.names=F, quote=""
)

samp_group_df <- read.table(
    argv[2], header=T, sep="\t",
    check.names=F, quote=""
)

colnames(samp_group_df) <- c("sample", "group")
SampNames <- samp_group_df$sample
ColNames <- samp_group_df$group

grp_color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")

source_color <- c("#8dd3c7", "#fdb462", "#b3de69", "#00688b", "#d9d9d9",
                  "#c34b48", "#fffac8", "#a163ad", "#ff1493", "#00c1ff",
                  "#a87c00", "#619cff", "#f8766d", "#008000", "#964500",
                  "#f6b3c3", "#5b8cb2", "#d2d28e", "#b297c9", "#ccebc5",
                  "#cc3a21", "#00ba38", "#a15d98", "#004a7f", "#ffc0cb",
                  "#848482", "#00ffff", "#ffed6f", "#9b9b9b", "#ffffe5",
                  "#bebada", "#fb8072", "#e6a141", "#80b1d3", "#f63200",
                  "#fccde5", "#bc80bd", "#5da5b3", "#ff00ff", "#556b2f",
                  "#008080", "#b41e00", "#6c7c32", "#cc8899", "#8f00ff",
                  "#ffa500", "#3b3b3b", "#ff0000", "#006400", "#00a3a3",
                  "#cd950c", "#00688b", "#8b795e", "#458b74", "#9acdf2",
                  "#ee9a00", "#00ff00", "#bebebe", "#9ac73f", "#beaed4")

ann_colors <- list(Group=generate_named_vector(unique(ColNames), grp_color))
annotation_col <- data.frame(Group=factor(ColNames))
rownames(annotation_col) <- SampNames

if('Source' %in% colnames(df)){
    Uniq_Source <- sort(unique(df$Source))
    annotation_row <- data.frame(Source = factor(df$Source))
    ann_colors$Source <- generate_named_vector(Uniq_Source,source_color)
    rownames(annotation_row) <- rownames(df)
}else{
    annotation_row<-NA
}

df <- subset(df,select=SampNames)
rownames(annotation_col)=names(df)

pdf(paste(argv[3], '/abundance_heatmap.pdf', sep=''), width=10, height=8)
pheatmap(
    log10(df+1),
    cluster_row=TRUE,
    cluster_col=FALSE,
    show_rownames=FALSE,
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
dev.off()
