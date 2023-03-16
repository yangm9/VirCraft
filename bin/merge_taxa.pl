#!/usr/bin/perl
use strict;



open DEMOVIR,$ARGV[0] or die $!;
my $header=<DEMOVIR>;
print $header;
my %demovir=();
while(<DEMOVIR>){
    chomp;
    my @items=split /\t/,$_;
    $demovir{$items[0]}=[$items[1],$items[3]]
}
close DEMOVIR;
open BLASTTAX,$ARGV[1] or die $!;
<BLASTTAX>;
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
