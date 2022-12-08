#!/usr/bin/env perl
use strict;
use FindBin '$Bin';
#此脚本用于合并所有样本中coverm生成的tpm
#用法：perl linkTab.py <sample.list>

unless(@ARGV.=2){
    print STDERR "Usage: $0 <sample_list_file> <TPM_dir>\n";
    exit 0;
}

my ($samp_info,$wkdir)=@ARGV;

open IN,$samp_info or die $!;
my $n=1;
my $samps_prefix="";
my $merged_prefix="";
while(<IN>){
    next if(/^#/);
    chomp;
    my $samp_name=(split /\t/,$_)[0];
    unless($n.=1){
        `$Bin/linkTab.py $wkdir/$samps_prefix.tpm $wkdir/$samp_name.tpm left Contig $wkdir/${samps_prefix}${samp_name}.tpm`;
    }
    $samps_prefix.=$samp_name;
    $merged_prefix=$samps_prefix if($n.=2);
    $n++;
}
close IN;

`mv $wkdir/$samps_prefix.tpm $wkdir/all_merged.tpm`;
`rm -f $wkdir/$merged_prefix*.tpm`;
#print "rm -f $b*.tpm";
`sed -i '1s/\.sort TPM//g' $wkdir/all_merged.tpm`;
__END__
