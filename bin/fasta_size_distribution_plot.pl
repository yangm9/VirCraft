#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# ---- 默认参数 ----
my $bin_width = 200;
my $log_y     = "FALSE";
my $help;

# ---- 解析命令行参数 ----
GetOptions(
    'bin-width=i' => \$bin_width,
    'log-y=s'     => \$log_y,
    'help|h'      => \$help
) or die "Error in command line arguments\n";

if ($help || @ARGV < 2) {
    print "Usage: perl $0 [--bin-width N] [--log-y TRUE|FALSE] <fasta1> [<fasta2> ...] <output_file>\n";
    exit;
}

# 最后一个参数是输出文件，其余是输入 fasta
my $output_pdf = pop @ARGV;
my @fasta_files = @ARGV;

# 临时中间文件
my $tmp_data = "temp_lengths_" . time() . ".csv";
my $tmp_r    = "temp_plot_" . time() . ".R";

# ==========================================
# 1. Perl 提取长度 (精准解析)
# ==========================================
open OUT, '>', $tmp_data or die "Cannot write to $tmp_data: $!";
print OUT "Sequence_ID,Length,Source\n";

$/ = ">";
foreach my $fasta(@fasta_files){
    if(!-e $fasta){
        warn "File $fasta not found, skipping...\n";
        next;
    }
    my $name = basename($fasta);
    open IN, $fasta or die "Cannot open $fasta: $!";
    my $current_len = 0;
    <IN>;
    while(<IN>){
        chomp;
        my($id, $seq) = split /\n/, $_, 2;
        $seq=~s/\s+//g; # 这里的正则确保只计算序列字符，排除换行、空格和数字
        my $length = length($seq); 
        print OUT "$id,$length,$name\n";
    }
    close IN;
}
$/ = "\n";
close OUT;

# ==========================================
# 2. 生成 R 脚本
# ==========================================
open RFILE, '>', $tmp_r or die "Cannot write to $tmp_r: $!";
print RFILE <<EOF;
library(ggplot2)

df <- read.csv("$tmp_data")
# 确保数据不为空
if (nrow(df) == 0) {
    stop("No sequence data extracted. Check your FASTA files.")
}

p <- ggplot(df, aes(x = Length, fill = Source)) +
    geom_histogram(binwidth = $bin_width, position = "stack", color = "white") +
    theme_bw() +
    labs(
        title = "Sequence Length Distribution",
        x = "Sequence length (bp)",
        y = "Sequence counts",
        fill = "Source file"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid.major = element_blank()
    )

# 处理对数轴逻辑
if (toupper("$log_y") == "TRUE") {
    p <- p + scale_y_log10()
}

ggsave("$output_pdf", p, width = 10, height = 6)
EOF
close RFILE;

# ==========================================
# 3. 执行 R 并清理
# ==========================================
system("Rscript $tmp_r");

# 删除中间文件
unlink($tmp_data, $tmp_r);

print "Successfully plotted distribution to: $output_pdf\n";
