#!/usr/bin/env perl
use strict;
use FindBin '$Bin';
#此脚本用于合并所有样本中coverm或salmon生成的tpm

unless(@ARGV>2){
    print STDERR "Usage: $0 <sample_list_file> <count_dir> <postfix> [Contig/Gene]\n";
    exit 0;
}
# Initialize the parameters 
my ($samp_info,$wkdir,$postfix,$object)=@ARGV;
$postfix||='tpm';
$object||='Contig';
my $topleft='Contig';
$topleft='Name' if($object eq 'Gene');
my %abd_hash=(
    'tpm'=>'TPM',
    'cov'=>'Mean'
);

# Define the variables
my $samp_number=(split /\s+/,`wc -l $samp_info`)[0];
my $samps_prefix="";
my $merged_prefix="";

# For Only one sample

if($samp_number==2){
    my $samp_name=(split /\t/,`sed -n '2p' $samp_info`)[0];
    `sed '1s/TPM/$samp_name/' $wkdir/${samp_name}_gene_quant/quant.sf|cut -f 1,4 > $wkdir/$samp_name.$postfix` if($object eq 'Gene');
    $samps_prefix=$merged_prefix=$samp_name;
    goto LABEL;
}

# For more than one sample

open IN,$samp_info or die $!;
my $n=1;
<IN>;
while(<IN>){
    chomp;
    my $samp_name=(split /\t/,$_)[0];
    `sed '1s/TPM/$samp_name/' $wkdir/${samp_name}_gene_quant/quant.sf|cut -f 1,4 > $wkdir/$samp_name.$postfix` if($object eq 'Gene');
    if($n==1){
        $samps_prefix=$samp_name;
        $merged_prefix=$samp_name;
    }else{
        #print "$Bin/linkTab.py $wkdir/$samps_prefix.$postfix $wkdir/$samp_name.$postfix left $topleft $wkdir/${samps_prefix}${n}.$postfix\n";
        `$Bin/linkTab.py $wkdir/$samps_prefix.$postfix $wkdir/$samp_name.$postfix left $topleft $wkdir/${samps_prefix}${n}.$postfix`;
        $samps_prefix.=$n;
        `rm -f $wkdir/$samp_name.$postfix`;
    }
    $n++;
}
close IN;

LABEL: `mv $wkdir/$samps_prefix.$postfix $wkdir/all_merged.$postfix`;
#`rm -f $wkdir/$merged_prefix*.$postfix`;
`sed -i '1s/\.sort $abd_hash{$postfix}//g' $wkdir/all_merged.$postfix` if($object eq 'Contig');
__END__
