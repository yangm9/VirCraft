#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my ($vrhyme_dir, $out_bin_file) = @ARGV;
die "Usage: $0 <vrhyme_dir> <output_fasta>\n" unless @ARGV == 2;

my @fasta_files = glob("$vrhyme_dir/vRhyme_best_bins_fasta/vRhyme_bin_*.fasta");
open OUT_BIN, ">$out_bin_file" or die $!;
foreach my $file(@fasta_files){
    my($joined_id, $joined_seq) = join_fasta_with_Ns($file);
    my $bin_name = basename($file);
    $bin_name =~ s/\.fasta$//i;
    print OUT_BIN ">$bin_name $joined_id\n$joined_seq\n";
}
close OUT_BIN;

sub join_fasta_with_Ns{
    my($fasta_file) = @_;
    local $/ = ">";
    open FASTA, $fasta_file or die "can not open file $fasta_file: $!";
    my(@ids, @seqs);
    <FASTA>;
    while(<FASTA>){
        chomp;
        my($id, $seq) = split /\n/, $_, 2;
        $id = (split /__/, $id, 2)[1];
        $seq =~ s/\s+//g;
        push @ids, $id;
        push @seqs, $seq;
    }
    close FASTA;
    local $/ = "\n";
    my $N_seq = 'N' x 1500;
    my $joined_id = join "|", @ids;
    my $joined_seq = join $N_seq, @seqs;
    return $joined_id, $joined_seq;
}
__END__
