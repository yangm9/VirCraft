#!/usr/bin/env perl
use strict;
use warnings;

# 输入和输出文件可通过参数传入，或使用默认值
my $input_file  = shift || 'genomad_output.tsv';
my $output_file = shift || 'genomad_taxa_parsed.tsv';

open(IN, '<', $input_file) or die "Cannot open input file $input_file: $!";
open(OUT, '>', $output_file) or die "Cannot open output file $output_file: $!";

# 打印表头
print OUT "Sequence_ID\tSuperrealm\tRealm\tKingdom\t\tClass\tOrder\tFamily\tGenus\tSpecies\n";

# 读取表头行
my $header = <IN>;
chomp $header;
my @columns = split /\t/, $header;

# 找到目标列索引
my %col_idx;
@col_idx{@columns} = (0..$#columns);

while(<IN>){
    chomp;
    my @fields = split /\t/;
    my $seq_id   = $fields[$col_idx{'seq_name'}];
    my $taxonomy = $fields[$col_idx{'taxonomy'}];

    # 处理 taxonomy 字段
    my ($superrealm, $realm, $kingdom, $phylum, $class, $order, $family, $genus, $species) = split /;/, $taxonomy;
    
    $superrealm ||= 'NA';
    $realm ||= 'NA';
    $kingdom ||= 'NA'; 
    $phylum ||= 'NA';
    $class ||= 'NA';
    $order ||= 'NA';
    $family ||= 'NA';
    $genus ||= 'NA';
    $species ||= 'NA';

    print OUT "$seq_id\t$realm\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
}

close IN;
close OUT;
__END__
