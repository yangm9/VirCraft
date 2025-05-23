#!/usr/bin/env Rscript
#old_name: stacked_multi_bar_plot.R
#yangm@idsse.ac.cn

argv <- commandArgs(T)

if(length(argv) < 2){
    stop("inputs: <viral_lifetype_quality_tsv> <barplot_dir>\n")
}

library(ggplot2)
library(ggnewscale)
library(RColorBrewer)

#类似于拉链功能，生成命名向量
generate_named_vector <- function(keys, values) {
    n <- length(keys)
    named_vector <- setNames(values[1:n], keys)
    return(named_vector)
}

# 数据框
dt <- read.table(
    argv[1], header = TRUE, sep = "\t",
    row.names = 1, check.names = F, quote = ""
)
dt$type[dt$type == "" | is.na(dt$type)] <- "Undetermined"

type_df <- data.frame(table(dt[, "type"]))
quality_df <- data.frame(table(dt[, "checkv_quality"]))

# 添加类别数量到类别名称中
type_df$LifeStyle <- paste(type_df$Var1, " (", type_df$Freq, ")", sep = "")
quality_df$Quality <- paste(quality_df$Var1, " (", quality_df$Freq, ")", sep = "")

colnames(type_df) <- c("Var1", "count", "LifeStyle")
colnames(quality_df) <- c("Var1", "count", "Quality")

# 计算百分比
type_df <- transform(type_df, percentage = count / sum(count) * 100)
quality_df <- transform(quality_df, percentage = count / sum(count) * 100)

# 从小到大排序
type_df <- type_df[order(type_df[, "percentage"]),]
quality_df <- quality_df[order(quality_df[, "percentage"]),]

# 颜色顺序
#sorted_Quality<-c("Complete","High-quality","Medium-quality","Low-quality","Not-determined")
sorted_Quality <- quality_df$Quality
#sorted_LifeStyle<-c("lysogenic","lytic","Undetermined")
sorted_LifeStyle <- type_df$LifeStyle

# 颜色设置
color_Quality <- generate_named_vector(sorted_Quality, brewer.pal(9, "Reds")[c(1, 3, 5, 7, 9)])
color_LifeStyle <- generate_named_vector(sorted_LifeStyle, brewer.pal(9, "Blues")[c(3, 6, 9)])

# 绘制堆叠条形图
barplot <- ggplot() +
    geom_bar(data = type_df, aes(x = "LifeStyle", y = percentage, fill = factor(LifeStyle, levels = sorted_LifeStyle)),
             stat = "identity", alpha = 0.7, position = "fill") +
    scale_fill_manual(values = color_LifeStyle,
                      breaks = sorted_LifeStyle,
                      name = "LifeStyle") +
    new_scale_fill() + # Define scales before initiating a new one
    geom_bar(data = quality_df, aes(x = "Quality", y = percentage, fill = factor(Quality, levels = sorted_Quality)),
             stat = "identity", alpha = 0.7, position = "fill") +
    scale_fill_manual(values = color_Quality,
                      breaks = sorted_Quality,
                      name = "Quality") +
    labs(x = "", y = "Percentage") +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black")
    )

pdf(paste(argv[2], '/phage_lifetype_quality_bar_plot.pdf', sep = ''), width = 5, height = 10)
barplot
dev.off()
