#color <- c(
#  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
#  "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
#  "#ccebc5", "#ffed6f", "#8C9EFF", "#9b9b9b", "#ffffe5",
#  "#fdb462", "#bebada", "#fb8072", "#80b1d3", "#ffffb3",
#  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"，
#  "#9ACD32", "#EECFA1", "#4F94CD", "#8B636C", "#EE9A00",
#  "#FF00FF", "#EECFA1", "#0000FF", "#BF3EFF", "#00FFFF",
#  "#A52A2A", "#FFFF00", "#FF00FF", "#00FF00", "#BEBEBE",
#  "#FF0000", "#EE9A00", "#006400", "#556B2F", "#556B2F",
#  "#CD950C", "#00688B", "#8B795E", "#458B74", "#9ACDF2"

#)

grp_color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#9edae5",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#17becf", "#dbdb8d", "#c7c7c7")

src_color <- c("#8dd3c7", "#fdb462", "#b3de69", "#00688b", "#d9d9d9",
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

# 创建一个1行n_cols列的矩阵，用于表示颜色条
n_cols <- length(color)
col_matrix <- matrix(1:n_cols, ncol = n_cols)

# 绘制颜色条
par(mar = c(5, 1, 1, 1))
image(col_matrix, col = color, axes = FALSE, xlab = "", ylab = "",
      main = "60 Colors", asp = 1)
