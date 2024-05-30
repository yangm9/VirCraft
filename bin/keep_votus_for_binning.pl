#!/usr/bin/perl
use strict;

die "$0 <checkv_quality_file> <vOTU_fasta> <fasta_for_binning>\n" if(@ARGV<3);

open QUAL, $ARGV[0] or die $!;
<QUAL>;
my %complete=();
while(<QUAL>){
    chomp;
    my @items=split /\t/,$_;
    if($items[7] eq 'Complete'){
        $complete{$items[0]}=1;
    }
}
close QUAL;

open VOTUS, $ARGV[1] or die $!;
open OUT, ">$ARGV[2]" or die $!;
$/='>';
<VOTUS>;
while(<VOTUS>){
    chomp;
    my($id,$seq)=split /\n/,$_,2;
    my @ids=split /\s+/,$id;
    if(@ids==1 and (!exists $complete{$id})){
        print OUT ">$id\n$seq";
    }
}
close VOTUS;
close OUT;
