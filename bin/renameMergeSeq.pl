#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;

=head1 renameMergeSeq.pl
Name:
    renameMergeSeq.pl
    Used Name: renameMAGsSeq.pl
Purpose: 
    1. Rename the sequence of MAGs/vOTUs to the format like "><File_name>_<seq_name>";
    2. Merge all sequence from the FastA files from a certain directory.
    Note: usually used for rename the MAGs and OTUs.
Usage:
    perl renameMergeSeq.pl -i <FastA_dir> [options] -o <merged.fa>
        -i <FastA_dir>: stored the FastA files with the postfix of ".fa"
        -o <merged.fa>: the output FastA file in which the sequence has been renamed
        -p <postfix>: the postfix of FastA file, i.e. fa, fasta, fna or faa, etc [default='fa']
        -t <seqtype>: the format of name, MAG or Seq, [default='Seq']
Author:
    Ming Yang yangm@idsse.ac.cn
Vwesion:
    0.0.2 2022.6.17 18:08
=cut

my ($fadir, $outfa, $postfix, $seqtype, $help);
GetOptions(
    "i:s" => \$fadir,
    "o:s" => \$outfa,
    "p:s" => \$postfix,
    "t:s" => \$seqtype,
    "help|h" => \$help
);
if($help || !defined $fadir || !defined $outfa){
    print STDERR `pod2text $0`;
    exit 0;
}

$postfix ||= "fa";
$seqtype ||= "Seq";
my @mags = glob "$fadir/*.$postfix";
open OUT, ">$outfa" or die $!;
foreach my $i(0 .. $#mags){
    my $mag_name = basename($mags[$i]);
    $mag_name =~ s/\.$postfix$//; #mag file name
    $mag_name =~ s/-//;
    open FA, $mags[$i] or die $!;
    while(<FA>){
        chomp;
        my $line = $_;
        if(/^>/){
            $line =~ s/^>//;
            $line = &fabricateMAGName($line) if($seqtype eq "MAG");
            $line =~ s/_length_.*_cov_.*$//;
            $line = $mag_name . "_" . $line; #For MAGs
            print OUT ">$line\n";
        }else{
            print OUT "$line\n";    
        }
    }
    close FA;
}
close OUT;

sub fabricateMAGName{
    my $line = shift;
    my @items = split /_/, $line;
    $line = join "_", @items[0, 1];
    return $line;
}
