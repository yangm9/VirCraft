#!/usr/bin/env Rscript

argv <- commandArgs(T)

if(length(argv) < 3){
    stop("inputs: <otus_tab.tsv> <group.txt> <outdir>")
}

if(!require("vegan")){install.packages("vegan")}
if(!require("ggplot2")){install.packages("ggplot2")}

library(vegan)   # 计算距离矩阵
library(ggplot2) # 绘图

otu <- read.table(
    argv[1], sep = "\t", header = TRUE,
    check.names = FALSE, row.names = 1
)

otu_t <- t(otu)
otu_dist <- vegdist(otu_t, method = "bray")

# 执行PCoA分析 ---------------------------
pcoa_res <- cmdscale(otu_dist, k = 2, eig = TRUE)

# 提取PCoA坐标 ---------------------------
df_points <- as.data.frame(pcoa_res$points)
colnames(df_points) <- c("PCoA1", "PCoA2")
df_points$sample <- row.names(df_points)

# 读取分组信息 ---------------------------
group <- read.table(argv[2], sep = "\t", header = TRUE)
group <- group[, 1:2]
colnames(group) <- c("sample", "group")

# 合并数据框 ---------------------------
df <- merge(df_points, group, by = "sample")

# 颜色方案 ---------------------------
color <- c("#1597A5", "#FFC24B", "#85B22E", "#FEB3AE", "#5F80B4", "#964500", "#619cff", "#8f00ff")

# 计算方差解释率 ---------------------------
eig_values <- pcoa_res$eig
var_explained <- eig_values / sum(eig_values) * 100

# 绘制PCoA图 ---------------------------
pcoa_plot <- ggplot(df, aes(x = PCoA1, y = PCoA2)) +
    theme_bw() +
    geom_point(aes(color = group), shape = 19, size = 3) +
    theme(panel.grid = element_blank()) +
    geom_vline(xintercept = 0, lty = "dashed", linewidth = 1, color = "grey50") +
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = "grey50") +
    geom_text(
        aes(label = sample, y = PCoA2 + 0.03, x = PCoA1 + 0.03, color = group),
        vjust = 0,
        size = 3.5,
        show.legend = FALSE
    ) +
    stat_ellipse(
        geom = "polygon",
        level = 0.95,
        linetype = 2,
        linewidth = 0.5,
        aes(fill = group),
        alpha = 0.2
    ) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    labs(
        x = paste0("PCoA1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PCoA2 (", round(var_explained[2], 1), "%)"),
        title = "PCoA based on Bray-Curtis distance"
    ) +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
    )

# 显示PCoA图 ---------------------------
pdf(paste(argv[3], '/PCoA.pdf', sep=''), width=8, height=5)
pcoa_plot
dev.off()

#ADONIS分析
group_list <- df$group
adonis_df <- adonis2(otu_t~group_list, method="bray", perm=999)
adonis_df <- cbind(ADONIS=row.names(adonis_df), adonis_df)
adonis_tab <- paste(argv[3], '/ADONIS_PCoA.tsv', sep='')#,col.names=TRUE,row.names=FALSE)
write.table(adonis_df, adonis_tab, sep='\t')

#ANOSIM分析
otu_t_dist <- vegdist(otu_t)
otu_t_ano <- with(group, anosim(otu_t_dist, group))
par(mar = c(5, 5, 4, 2))
pdf(paste(argv[3], '/ANOSIM_PCoA.pdf', sep=''), width=10, height=8)
plot(otu_t_ano, xlab = "", ylab = "Dissimilarity Rank Value")
dev.off()
