#!/usr/bin/perl -w
use strict;

die "perl $0 <fasta_file> <cutted_fasta_prefix> <cutoff_len>" unless (@ARGV==3);
my $fa=$ARGV[0];
my $out=$ARGV[1];
my $cutoff=$ARGV[2];
open OUT,">$out.gt$cutoff.fa" || die "Can not open $out.fa\n";
open FA,$fa||die "Can not read $fa\n";
$/=">";
<FA>;
my $num=1;
while(<FA>){
    chomp;
    my($id,$seq)=split /\n/,$_,2;
    $id=(split /\s+/,$id)[0];
    $seq=~s/\n//g;
    if(length($seq)>=$cutoff) {
        $id=">".$id."_".$num;
        $seq=~s/(.*)/\U$1/; #upper
        print OUT $id."\n".$seq."\n";
        $num++;
    }
}
$/="\n";
close FA;
close OUT;
__END__
