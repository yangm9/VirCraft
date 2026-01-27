#!/usr/bin/env Rscript

# ==============================
# Stacked Length Distribution from FASTA
# Command-line version (no Biostrings)
# ==============================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ==============================
# 1. 解析命令行参数
# ==============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript plot_fasta_lengths.R ",
    "[--bin-width N] [--log-y TRUE|FALSE] ",
    "<fasta1> [<fasta2> ...] <output_file>"
  )
}

# ---- 默认参数 ----
bin_width <- 200
use_log_y <- FALSE

# ---- 解析可选参数 ----
parse_option <- function(flag, default) {
  idx <- match(flag, args)
  if (!is.na(idx)) {
    value <- args[idx + 1]
    args <<- args[-c(idx, idx + 1)]
    return(value)
  }
  default
}

bin_width <- as.numeric(parse_option("--bin-width", bin_width))
use_log_y <- as.logical(parse_option("--log-y", use_log_y))

if (is.na(bin_width) || bin_width <= 0) {
  stop("--bin-width must be a positive number")
}

if (is.na(use_log_y)) {
  stop("--log-y must be TRUE or FALSE")
}

# ---- 剩余参数 ----
output_file <- args[length(args)]
fasta_files <- args[-length(args)]

# ==============================
# 2. FASTA 读取函数（base R）
# ==============================
read_fasta_lengths <- function(fasta_file) {
  lines <- readLines(fasta_file, warn = FALSE)
  lines <- lines[nzchar(lines)]

  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", fasta_file)
  }

  start_idx <- header_idx + 1
  end_idx <- c(header_idx[-1] - 1, length(lines))

  seq_lengths <- mapply(
    function(start, end) {
      nchar(paste0(lines[start:end], collapse = ""))
    },
    start_idx,
    end_idx
  )

  data.frame(
    Length = seq_lengths,
    Source = basename(fasta_file),
    stringsAsFactors = FALSE
  )
}

# ==============================
# 3. 读取所有 FASTA
# ==============================
df <- do.call(
  rbind,
  lapply(fasta_files, read_fasta_lengths)
)

# ==============================
# 4. 绘图
# ==============================
p <- ggplot(df, aes(x = Length, fill = Source)) +
  geom_histogram(
    binwidth = bin_width,
    position = "stack",
    color = "white"
  ) +
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence length (bp)",
    y = "Sequence counts",
    fill = "FASTA file"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq(0, max(df$Length, na.rm = TRUE), by = 2000)
  )

if (use_log_y) {
  p <- p + scale_y_log10()
}

# ==============================
# 5. 保存图像
# ==============================
ggsave(output_file, p, width = 10, height = 6)
message("Plot saved to ", output_file)
