#!/usr/bin/env Rscript
#old_name: 
#yangm@idsse.ac.cn

argv <- commandArgs(T)

if(length(argv) < 2){
    stop("inputs: <quality_summary.tsv> <boxplot_dir>\n")
}

if(!require('ggplot2')){
    install.packages('ggplot2')
    library(ggplot2)
}
if(!require('ggpubr')){
    install.packages('ggpubr')
    library(ggpubr)
}
if(!require('car')){
    install.packages('car')
    library(car)
}

data <- read.table(argv[1], header=T, sep='\t')
all_contig_length <- data$contig_length

# 将 checkv_quality 转换为因子，并设置因子的级别顺序
data$checkv_quality <- factor(data$checkv_quality, 
                                levels = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))

# 获取最大值
ymax <- max(data$contig_length, na.rm = TRUE)

### 绘图
boxplot <- ggplot(data = data1, aes(x = checkv_quality, y = contig_length)) + 
    geom_jitter(aes(colour = checkv_quality), alpha = 0.3, size = 3) +
    geom_boxplot(alpha = 0.5, size = 0.8, width = 0.6, colour = "black") + 
    scale_color_manual(
        limits = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"), 
        values = c("#e6a141", "#85B22E", "#5F80B4", "#f8766d", "#d9d9d9")) + 
    stat_compare_means(
        comparisons = combn(levels(data1$checkv_quality), 2, simplify = FALSE),
        method = "wilcox.test",
        label = "p.signif",
        tip.length = 0.01,
        step.increase = 0.08,
        size = 4
    ) +
    theme_classic(base_line_size = 0.5) +
    labs(title = "", x = "", y = "Genome size") + 
    scale_y_continuous(limits = c(0, ymax * 2)) +
    theme(
        plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    )

pdf(paste(argv[2],'/quality_length_box_plot.pdf',sep=''),width=5,height=10)
boxplot
dev.off()
