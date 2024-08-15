#!/usr/bin/env Rscript
#coding:utf8
#yangm@idsse.ac.cn

# 加载所需的R包
if(!require('tools')) {install.packages('tools'); library(tools)}
if(!require('vegan')) {install.packages('vegan'); library(vegan)}
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('patchwork')) {install.packages('patchwork'); library(patchwork)}
if(!require('ggpubr')) {install.packages('ggpubr'); library(ggpubr)}

# 接收命令行位置参数
argv<-commandArgs(T)

# 判断命令行参数个数
if (length(argv)<2){
    stop("inputs: <all_abundance.xls> <samp_info.list> <alpha_diversity.xls>")
}

# 读取数据
otu <- read.table(argv[1], sep="\t", header=T,
                  check.names=FALSE, row.names=1)
otu_t <- t(otu)

# 计算多样性指数
Shannon = diversity(otu_t, "shannon") 
Simpson = diversity(otu_t, "simpson") 
InSimpson = diversity(otu_t, "inv")
specN = specnumber(otu_t)
Pielou = Shannon / log(specN)  # 计算Pielou均匀度指数

# 合并数据并导出结果
alpha_df = t(rbind(specN, Pielou, Shannon, Simpson, InSimpson))
alpha_df <- as.data.frame(alpha_df)
alpha_df$Sample <- rownames(alpha_df)
rownames(alpha_df) <- NULL
alpha_df <- alpha_df[, c(ncol(alpha_df), 1:(ncol(alpha_df) - 1))] # 将最后一列Sample移动到第一列
outprefix=file_path_sans_ext(basename(argv[1]))
alpha_tab = paste(argv[3],'/',outprefix,'.alpha_diversity.xls',sep='')
write.table(alpha_df, alpha_tab, sep='\t', row.names=FALSE)

# 读取样本分组信息并合并
samp_group_df <- read.table(argv[2], header=T, sep="\t", # samp.info.xls
                            check.names=F, quote="")
samp_group_df <- samp_group_df[, 1:2]
alpha_df <- merge(samp_group_df, alpha_df, by.x="Sample", by.y="Sample")

# 设置颜色
color1 <- c("#ff7f0e", "#2ca02c", "#d62728", "#9467bd","#5F80B4")
group_number <- length(unique(samp_group_df$Group))
color1 <- color1[1:group_number]

# 仅保留显著性结果的比较
get_significant_comparisons <- function(df, y_value) {
    comparisons <- combn(unique(df$Group), 2, simplify = FALSE)
    p_values <- sapply(comparisons, function(comp) {
        group1 <- df[df$Group == comp[1], y_value]
        group2 <- df[df$Group == comp[2], y_value]
        wilcox.test(group1, group2)$p.value
    })
    significant_comparisons <- comparisons[p_values < 0.05]
    return(significant_comparisons)
}

# 定义一个绘图函数，并添加显著性标记
plot_with_significance <- function(df, y_value, y_label) {
    significant_comparisons <- get_significant_comparisons(df, deparse(substitute(y_value)))
    ggplot(df, aes(x=Group, y={{ y_value }})) +
        geom_boxplot(outlier.size=1, fill=color1) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              panel.grid=element_blank(),
              panel.background=element_rect(fill='transparent',color='black')) +
        labs(x='', y = y_label) +
        stat_compare_means(method="wilcox.test", label="p.signif", 
                           comparisons=significant_comparisons, 
                           hide.ns=TRUE)  # 仅显示显著性结果的连线和标记
}

p1 <- plot_with_significance(alpha_df, specN, 'specN')
p2 <- plot_with_significance(alpha_df, Shannon, 'Shannon')
p3 <- plot_with_significance(alpha_df, Simpson, '1-Simpson')
p4 <- plot_with_significance(alpha_df, Pielou, 'Pielou')

# 保存结果
pdf(paste(argv[3],'/',outprefix,'.alpha_diversity.boxplot.pdf',sep=''),width=10,height=8)
p1 + p2 + p3 + p4
dev.off()
