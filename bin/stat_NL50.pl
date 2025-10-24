#!/usr/bin/perl -w
## zhushilin@genomics.org.sn 2011-05-20
## yangm@idsse.ac.cn 2023-11-25 

use strict;
die "$0 <scaffold.fa> <output> <cutoff>\n" if(@ARGV < 2);

my ($fasta, $output, $cutoff) = @ARGV;
$cutoff = 100 if (!defined $cutoff);

my (@scaffold_length, @contig_length);
my ($total_scaffold_length, $total_contig_length) = (0, 0);	### total_len
my ($scaffold_gc, $contig_gc) = (0, 0);
my ($scaf_gt_100, $cont_gt_100) = (0, 0);		### number (>=cutoff)
my ($scaf_gt_500, $cont_gt_500) = (0, 0);		### number (>=500)
my ($scaf_gt_1000, $cont_gt_1000) = (0, 0);		### number (>=1000)
my ($scaf_gt_1500, $cont_gt_1500) = (0, 0);		### number (>=1500)
my ($scaf_gt_2000, $cont_gt_2000) = (0, 0);		### number (>=2000)
my ($scaf_gt_5000, $cont_gt_5000) = (0, 0);		### number (>=5000)
my ($scaf_gt_10000, $cont_gt_10000) = (0, 0);		### number (>=10000)

open FILE, $fasta or die "Can not open the input file $fasta: $!";
$/ = ">";
<FILE>;
while(<FILE>){
	chomp;
	my($tag, $seq) = split /\n/, $_, 2;
	($tag) = $tag =~/^(\S+)/;
	$seq =~s/\W+//g;
	$seq = uc $seq;
	### scaffold
	my $scaf_len = length $seq;
	if ($scaf_len >= $cutoff){
		push @scaffold_length, $scaf_len;
		$total_scaffold_length += $scaf_len;
		my $this_scaf_gc += $seq =~tr/GC/GC/;
		$scaffold_gc += $this_scaf_gc;
		$scaf_gt_100++ if $scaf_len >= $cutoff;
		$scaf_gt_500++ if $scaf_len >= 500;
		$scaf_gt_1000++ if $scaf_len >= 1000;
		$scaf_gt_1500++ if $scaf_len >= 1500;
		$scaf_gt_2000++ if $scaf_len >= 2000;
		$scaf_gt_5000++ if $scaf_len >= 5000;
		$scaf_gt_10000++ if $scaf_len >= 10000;
	}
	### contig
	my @contig = split /N+/i, $seq;
	foreach my $cont_seq (@contig){
		my $cont_len = length $cont_seq;
		if ($cont_len >= $cutoff){
			push @contig_length, $cont_len;
			$total_contig_length += $cont_len;
			my $this_cont_gc += $cont_seq =~tr/GC/GC/;
			$contig_gc += $this_cont_gc;
			$cont_gt_100++ if $cont_len >= $cutoff;
			$cont_gt_500++ if $cont_len >= 500;
		    $cont_gt_1000++ if $cont_len >= 1000;
		    $cont_gt_1500++ if $cont_len >= 1500;
			$cont_gt_2000++ if $cont_len >= 2000;
			$cont_gt_5000++ if $cont_len >= 5000;
			$cont_gt_10000++ if $cont_len >= 10000;
		}
	}
}
$/ = "\n";
close FILE;

### deal scaffold
@scaffold_length = sort { $b <=> $a } @scaffold_length;
my $scaffold_max_length = $scaffold_length[0];		### max_len
my $scaffold_gc_rate = sprintf "%.3f", $scaffold_gc / $total_scaffold_length;	### GC_rate
my %scaffold_Nx_stats = calculate_Nx_stats(\@scaffold_length, $total_scaffold_length);	### statistic

### deal contig
@contig_length = sort { $b <=> $a } @contig_length;
my $contig_max_length = $contig_length[0];		### max_len
my $contig_gc_rate = sprintf "%.3f", $contig_gc/$total_contig_length;	### GC_rate
my %contig_Nx_stats = calculate_Nx_stats(\@contig_length, $total_contig_length);	### statistic

### print result
open OUT, ">$output" or die "Can not open the output file $output: $!";
print OUT "\tscaffold\t\tcontig\n";
print OUT "\tlength(bp)\tnumber\tlength(bp)\tnumber\n";
print OUT "max_len\t$scaffold_max_length\t-\t$contig_max_length\t-\n";
for(my $i=1; $i<=9; $i++){
	print OUT "N/L$i" . "0\t$scaffold_Nx_stats{$i}{l}\t$scaffold_Nx_stats{$i}{n}\t$contig_Nx_stats{$i}{l}\t$contig_Nx_stats{$i}{n}\n";
}
print OUT "Total_length\t$total_scaffold_length\t-\t$total_contig_length\t-\n";
print OUT "number>=$cutoff" . "bp\t-\t$scaf_gt_100\t-\t$cont_gt_100\n";
print OUT "number>=500bp\t-\t$scaf_gt_500\t-\t$cont_gt_500\n";
print OUT "number>=1000bp\t-\t$scaf_gt_1000\t-\t$cont_gt_1000\n";
print OUT "number>=1500bp\t-\t$scaf_gt_1500\t-\t$cont_gt_1500\n";
print OUT "number>=2000bp\t-\t$scaf_gt_2000\t-\t$cont_gt_2000\n";
print OUT "number>=5000bp\t-\t$scaf_gt_5000\t-\t$cont_gt_5000\n";
print OUT "number>=10000bp\t-\t$scaf_gt_10000\t-\t$cont_gt_10000\n";
print OUT "GC_rate\t$scaffold_gc_rate\t-\t$contig_gc_rate\t-\n";
close OUT;

#===================================================================
sub calculate_Nx_stats{
	my ($lengths_array_ref, $total_length) = @_;
	my %nx_statistics;
	my $cumulative_count = 0;
	my $cumulative_length = 0;
	foreach (@$lengths_array_ref){
		$cumulative_count++;
		$cumulative_length += $_;
		if(!exists $nx_statistics{1} && $cumulative_length >= $total_length * 0.1){
			$nx_statistics{1}{'l'} = $_;
			$nx_statistics{1}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{2} && $cumulative_length >= $total_length * 0.2){
			$nx_statistics{2}{'l'} = $_;
			$nx_statistics{2}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{3} && $cumulative_length >= $total_length * 0.3){
			$nx_statistics{3}{'l'} = $_;
			$nx_statistics{3}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{4} && $cumulative_length >= $total_length * 0.4){
			$nx_statistics{4}{'l'} = $_;
			$nx_statistics{4}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{5} && $cumulative_length >= $total_length * 0.5){
			$nx_statistics{5}{'l'} = $_;
			$nx_statistics{5}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{6} && $cumulative_length >= $total_length * 0.6){
			$nx_statistics{6}{'l'} = $_;
			$nx_statistics{6}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{7} && $cumulative_length >= $total_length * 0.7){
			$nx_statistics{7}{'l'} = $_;
			$nx_statistics{7}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{8} && $cumulative_length >= $total_length * 0.8){
			$nx_statistics{8}{'l'} = $_;
			$nx_statistics{8}{'n'} = $cumulative_count;
		}
		if(!exists $nx_statistics{9} && $cumulative_length >= $total_length * 0.9){
			$nx_statistics{9}{'l'} = $_;
			$nx_statistics{9}{'n'} = $cumulative_count;
			last;
		}
	}
	return %nx_statistics;
}
