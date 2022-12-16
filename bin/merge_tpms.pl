#!/usr/bin/env perl
use strict;
use FindBin '$Bin';
#此脚本用于合并所有样本中coverm生成的tpm
#用法：perl linkTab.py <sample.list>

unless(@ARGV==3){
    print STDERR "Usage: $0 <sample_list_file> <count_dir> <postfix>\n";
    exit 0;
}

my ($samp_info,$wkdir,$postfix)=@ARGV;
$postfix||='count';

open IN,$samp_info or die $!;
my $n=1;
my $samps_prefix="";
my $merged_prefix="";
while(<IN>){
    next if(/^Sample/);
    chomp;
    my $samp_name=(split /\t/,$_)[0];
    unless($n==1){
        `$Bin/linkTab.py $wkdir/$samps_prefix.$postfix $wkdir/$samp_name.$postfix left Contig $wkdir/${samps_prefix}${samp_name}.$postfix`;
    }
    $samps_prefix.=$samp_name;
    $merged_prefix=$samps_prefix if($n==2);
    $n++;
}
close IN;

`mv $wkdir/$samps_prefix.$postfix $wkdir/all_merged.$postfix`;
`rm -f $wkdir/$merged_prefix*.$postfix`;
#print "rm -f $b*.$postfix";
`sed -i '1s/\.sort TPM//g' $wkdir/all_merged.$postfix`;
__END__
