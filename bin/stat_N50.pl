#!/usr/bin/perl -w
## zhushilin@genomics.org.sn 2011-05-20
## yangm@idsse.ac.cn 2023-11-25 

use strict;
die "$0 <scaffold.fa> <output> <cutoff>\n" if(@ARGV <2);

my ($fasta,$output,$cutoff)=@ARGV;
$cutoff=100 if (!defined $cutoff);

my (@sca_len,@con_len);
my ($total_slen,$total_clen)=(0,0);	### total_len
my ($sca_gc,$con_gc)=(0,0);
my ($sca100,$con100)=(0,0);		### number (>=cutoff)
my ($sca500,$con500)=(0,0);		### number (>=500)
my ($sca1000,$con1000)=(0,0);		### number (>=1000)
my ($sca1500,$con1500)=(0,0);		### number (>=1500)
my ($sca2000,$con2000)=(0,0);		### number (>=2000)
my ($sca5000,$con5000)=(0,0);		### number (>=5000)
my ($sca10000,$con10000)=(0,0);		### number (>=10000)

open FILE,"<","$fasta";
$/=">";
<FILE>;
while (<FILE>){
	chomp;
	my ($tag,$seq)=split(/\n/,$_,2);
	($tag)=$tag=~/^(\S+)/;
	$seq=~s/\W+//g;
	$seq= uc $seq;
	### scaffold
	my $sl=length $seq;
	if ($sl>=$cutoff){
		push @sca_len,$sl;
		$total_slen+=$sl;
		my $this_s_gc+=$seq=~tr/GC/GC/;
		$sca_gc+=$this_s_gc;
		$sca100++ if $sl>=$cutoff;
		$sca500++ if $sl>=500;
		$sca1000++ if $sl>=1000;
		$sca1500++ if $sl>=1500;
		$sca2000++ if $sl>=2000;
		$sca5000++ if $sl>=5000;
		$sca10000++ if $sl>=10000;
	}
	### contig
	my @contig=split /N+/i,$seq;
	foreach my $one_c (@contig){
		my $cl=length $one_c;
		if ($cl>=$cutoff){
			push @con_len,$cl;
			$total_clen+=$cl;
			my $this_c_gc+=$one_c=~tr/GC/GC/;
			$con_gc+=$this_c_gc;
			$con100++ if $cl>=$cutoff;
			$con500++ if $cl>=500;
		    $con1000++ if $sl>=1000;
		    $con1500++ if $sl>=1500;
			$con2000++ if $cl>=2000;
			$con5000++ if $cl>=5000;
			$con10000++ if $cl>=10000;
		}
	}
}
$/="\n";
close FILE;

### deal scaffold
@sca_len=sort { $b <=> $a } @sca_len;
my $sca_num=scalar @sca_len;		### number
my $smax_len=$sca_len[0];		### max_len
my $sca_gc_rate=sprintf "%.3f",$sca_gc/$total_slen;	### GC_rate
my %sN=stalen(\@sca_len,$total_slen);	### statistic

### deal contig
@con_len=sort { $b <=> $a } @con_len;
my $con_num=scalar @con_len;		### number
my $cmax_len=$con_len[0];		### max_len
my $con_gc_rate=sprintf "%.3f",$con_gc/$total_clen;	### GC_rate
my %cN=stalen(\@con_len,$total_clen);	### statistic
### print result
open OUT,">$output" || die "Fail to open output file: N50.xls\n";
print OUT "\tscaffold\t\tcontig\n";
print OUT "\tlength(bp)\tnumber\tlength(bp)\tnumber\n";
print OUT "max_len\t$smax_len\t-\t$cmax_len\t-\n";
for(my $i=1;$i<=9;$i++){
	print OUT "N$i"."0\t$sN{$i}{l}\t$sN{$i}{n}\t$cN{$i}{l}\t$cN{$i}{n}\n";
}
print OUT "Total_length\t$total_slen\t-\t$total_clen\t-\n";
print OUT "number>=$cutoff"."bp\t-\t$sca100\t-\t$con100\n";
print OUT "number>=500bp\t-\t$sca500\t-\t$con500\n";
print OUT "number>=1000bp\t-\t$sca1000\t-\t$con1000\n";
print OUT "number>=1500bp\t-\t$sca1000\t-\t$con1500\n";
print OUT "number>=2000bp\t-\t$sca2000\t-\t$con2000\n";
print OUT "number>=5000bp\t-\t$sca5000\t-\t$con5000\n";
print OUT "number>=10000bp\t-\t$sca10000\t-\t$con10000\n";
print OUT "GC_rate\t$sca_gc_rate\t-\t$con_gc_rate\t-\n";
close OUT;

#===================================================================
sub stalen{
	my ($this_len,$this_total)=@_;
	my %this_hash;
	my $num_temp=0;
	my $len_temp=0;
	foreach (@$this_len){
		$num_temp++;
		$len_temp+=$_;
		if (!exists $this_hash{1} && $len_temp>=$this_total*0.1){
			$this_hash{1}{'l'}=$_;
			$this_hash{1}{'n'}=$num_temp;
		}
		if (!exists $this_hash{2} && $len_temp>=$this_total*0.2){
			$this_hash{2}{'l'}=$_;
			$this_hash{2}{'n'}=$num_temp;
		}
		if (!exists $this_hash{3} && $len_temp>=$this_total*0.3){
			$this_hash{3}{'l'}=$_;
			$this_hash{3}{'n'}=$num_temp;
		}
		if (!exists $this_hash{4} && $len_temp>=$this_total*0.4){
			$this_hash{4}{'l'}=$_;
			$this_hash{4}{'n'}=$num_temp;
		}
		if (!exists $this_hash{5} && $len_temp>=$this_total*0.5){
			$this_hash{5}{'l'}=$_;
			$this_hash{5}{'n'}=$num_temp;
		}
		if (!exists $this_hash{6} && $len_temp>=$this_total*0.6){
			$this_hash{6}{'l'}=$_;
			$this_hash{6}{'n'}=$num_temp;
		}
		if (!exists $this_hash{7} && $len_temp>=$this_total*0.7){
			$this_hash{7}{'l'}=$_;
			$this_hash{7}{'n'}=$num_temp;
		}
		if (!exists $this_hash{8} && $len_temp>=$this_total*0.8){
			$this_hash{8}{'l'}=$_;
			$this_hash{8}{'n'}=$num_temp;
		}
		if (!exists $this_hash{9} && $len_temp>=$this_total*0.9){
			$this_hash{9}{'l'}=$_;
			$this_hash{9}{'n'}=$num_temp;
			last;
		}
	}
	return %this_hash;
}
