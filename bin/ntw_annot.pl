#!/usr/bin/env perl

#NTW为vContact2生成的蛋白质相互关系的结果。
#本程序的功能：过滤掉NTW中与研究无关的行并添加环境信息。

use strict;
unless(@ARGV==4){
    print STDERR "Usage: $0 <ntw> <ntw_anno_file> <Seq_Name_Token> <Currnet Environment>\n"
}

open NTW,$ARGV[0] or die $!;
open ANNO,">$ARGV[1]" or die $!;
print ANNO "Name\tTaget\tRelationship\tEnvironment\n";
while(<NTW>){
    chomp;
    next unless(/^$ARGV[2]/);
    next if(/nanopore/); #该行需要注释掉
    my @items=split /\s+/,$_;
    my $env="";
    if($items[1]=~/~/){
        $env="vContact2 Reference";
    }elsif($items[1]=~/^$ARGV[2]/){
        $env=$ARGV[3];
    }elsif($items[1]=~/^TS/){
        $env="Deep Sea (Depth>5000m)";
    }else{
        $env=(split /_/,$items[1])[0];
    }
    push @items,$env;
    my $line=join "\t",@items;
    print ANNO "$line\n";
}
close NTW;
close ANNO;
