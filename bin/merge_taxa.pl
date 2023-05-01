#!/usr/bin/perl
use strict;

=head1 merge_taxa.pl
Name:
    merge_taxa.pl
Purpose:
    Merge the classification results of NCBI BLAST and DemoVir.
Usage:
    perl merge_taxa.pl <DemoVir_Taxa> <NCBI_BLAST_Taxa> > <Merged_Taxa>
Author:
    Ming Yang yangm@idsse.ac.cn
Vwesion:
    0.0.1 2023.3.16 16:28
=cut

unless(@ARGV==2){
    print STDERR `pod2text $0`;
    exit 0;
}

open DEMOVIR,$ARGV[0] or die $!;
my $header=<DEMOVIR>;
my %demovir=();
while(<DEMOVIR>){
    chomp;
    my @items=split /\t/,$_;
    $demovir{$items[0]}=[$items[1],$items[3]]
}
close DEMOVIR;
open BLASTTAX,$ARGV[1] or die $!;
my $header=<BLASTTAX>;
print $header;
while(<BLASTTAX>){
    chomp;
    my @items=split /\t/,$_;
    print "$items[0]\t";
    if($items[1] eq 'Unassigned'){
        if(exists $demovir{$items[0]}){
            print "$demovir{$items[0]}[0]\t";
        }else{
            print "$items[1]\t";
        }
    }else{
        print "$items[1]\t";    
    }
    if($items[2] eq 'Unassigned'){
        if(exists $demovir{$items[0]}){
            print "$demovir{$items[0]}[1]\n";
        }else{
            print "$items[2]\n";
        }
    }else{
        print "$items[2]\n";    
    }
}
close BLASTTAX;
__END__
