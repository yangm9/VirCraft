#!/usr/bin/perl

=head1 Description
    Purpose: Remove contigs with a length < 5,000 bp, identified as a prophage, or having a CheckV quality of 'Complete'.
             Additionally, output filtered (removed) contigs (CheckV Complete or <5000bp or prophage=Yes)
    Useage: $0 <checkv_quality_file> <viral_contig_fasta> <fasta_for_binning>
    Auther: Ming Yang(yangming@genomics.cn)
    Version: 0.0.1 2023-3-21 21:45
=cut

use strict;
use Term::ANSIColor;

if(@ARGV < 3){
    print STDERR color("yellow"), `pod2text $0`, color("reset");    
    exit 0;
}

# Mark complete viral contigs, proviral contigs, and binning-ready contigs
my(%complete, %provirus, %vctg_for_binning);
open QUAL, $ARGV[0] or die $!;
<QUAL>;
while(<QUAL>){
    chomp;
    my @items = split /\t/, $_;
    if($items[7] eq 'Complete'){
        $complete{$items[0]} = 1; # complete viral contigs
    }elsif($items[2] eq 'Yes'){
        $provirus{$items[0]} = 1; #proviral contigs
    }elsif($items[1] >= 5000){
        $vctg_for_binning{$items[0]} = 1; #binning-ready contigs
    }
}
close QUAL;

# Output complete viral contigs, proviral contigs, and binning-ready contigs
open VCTGS, $ARGV[1] or die $!;
open COMPLETE, ">$ARGV[2]_vctg_complete.fa" or die $!;
open PROVIRUS, ">$ARGV[2]_vctg_provirus.fa" or die $!;
open OUT, ">$ARGV[2]_vctg_for_binning.fa" or die $!;
$/ = '>';
<VCTGS>;
while(<VCTGS>){
    chomp;
    my($id, $seq) = split /\n/, $_, 2;
    my @ids = split /\s+/, $id;
    if(exists $complete{$ids[0]}){
        print COMPLETE ">$id\n$seq";
    }elsif(exists $provirus{$ids[0]}){
        print PROVIRUS ">$id\n$seq";
    }elsif(exists $vctg_for_binning{$ids[0]}){
        print OUT ">$id\n$seq";
    }else{
        print "$id\n";
    }
}

close VCTGS;
close COMPLETE;
close PROVIRUS;
close OUT;
