#!/usr/bin/perl -w

use strict;

die <<"USAGE" unless (@ARGV == 3 || @ARGV == 4);
Usage:
  perl $0 <fasta_file> <output_prefix> <min_len_len>
  perl $0 <fasta_file> <output_prefix> <min_len> <max_len>

Examples:
  perl $0 input.fa output 300
  perl $0 input.fa output 300 1000
  perl $0 input.fa output 0 500
USAGE

my($fa, $out, $min_len, $max_len);

if(@ARGV == 3){
    ($fa, $out, $min_len) = @ARGV;
    $max_len = 0;  # Disable the upper limit
}else{
    ($fa, $out, $min_len, $max_len) = @ARGV;
}

my $outfile;
if($min_len == 0){
    $outfile = "$out.lt$max_len.fa";
}elsif($max_len) {
    $outfile = "$out.${min_len}_$max_len.fa";
}else{
    $outfile = "$out.gt$min_len.fa";
}

open FA, $fa || die "Can not read $fa\n";
open OUT, ">$outfile" || die "Can not open $outfile\n";

$/ = ">";
<FA>;
my $num = 1;
while(<FA>){
    chomp;
    my($id,$seq) = split /\n/, $_, 2;
    $id = (split /\s+/, $id)[0];
    $seq =~ s/\n//g;
    my $len = length($seq);
    if($min_len == 0){
        next unless $len <= $max_len;
    }elsif($max_len){
        next unless $len >= $min_len && $len <= $max_len;
    }else{
        next unless $len >= $min_len;
        $id .= "_$num";
    }
    $seq =~ s/(.*)/\U$1/; # upper
    print OUT ">$id\n$seq\n";
    $num++;
}
$/ = "\n";

close FA;
close OUT;
__END__
