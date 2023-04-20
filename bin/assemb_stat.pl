#!/usr/bin/perl
use strict;
use Getopt::Long;
my ($stat_lv,$clenf,$getcon,$item,$slenf,$cuta,$cutb,$nnum,$help);
GetOptions(
        "a:i"=>\$cuta,
        "b:i"=>\$cutb,
        "l:s"=>\$stat_lv,
        "c:s"=>\$clenf,
        "s:s"=>\$slenf,
        "n"=>\$nnum,
        "m:s"=>\$item,
        "g"=>\$getcon,
        "h"=>\$help
);
if(@ARGV !=2 || $help){
	die"Name: assemb_stat.pl
Describe: to statistics the assemble sequence length information
Author: liuwenbin, liuwenbin\@genomic.org.cn
Version: 1.0, Date: 2011-03-07
Usage: perl assemb_stat.pl <cont.fa> <scaf.fa>  >stat.tab
    -a <n>       conting length cout off, default=100
    -b <n>       scaffold length cout off, default=100
    -l <s>       Nxx% statistic leval, default=90,80,70,60,50
    -c <s>       conting sequence len file, default caculated by process
    -s <s>       scaffold sequence len file, default caculated by process
    -m <s>       item name, default=Contig,Scaffold
    -g           get contig len from id line
    -n           out put N_num at len out file
    -h           output help information to screen
Note: when -c or -s infile not exist, len result will output to it\n";
}
my ($cont,$scaf) = @ARGV;
$stat_lv ||= "90,80,70,60,50";
$getcon ||= 0;
my @slv = reverse(sort {$a<=>$b}(split/,/,$stat_lv));
my @sign;
foreach(@slv){push @sign,"N$_";}
push @sign,('Longest','Toal_Size','Number(>100bp)','Number(>2kb)');
my (@size1,@num1,@size2,@num2);
my @items = $item ? split/,/,$item : qw(Contig Scaffold);
single_stat($cont,\@slv,\@size1,\@num1,$clenf,$cuta,$getcon);#sub3
single_stat($scaf,\@slv,\@size2,\@num2,$slenf,$cutb,0,$nnum);#sub3
print "\t\t\t$items[0]\t\t\t$items[1]\n\t\tSize(bp)\tNumber\tSize(bp)\tNumber\n";
foreach my $i(0..$#sign){
	($i  < $#sign - 2) ? (print "$sign[$i]\t\t$size1[$i]\t\t$num1[$i]\t$size2[$i]\t\t$num2[$i]\n") :
	($i == $#sign - 2) ? (print "$sign[$i]\t$size1[$i]\t$num1[$i]\t$size2[$i]\t$num2[$i]\n") :
	(print "$sign[$i]\t$size1[$i]\t\t$num1[$i]\t$size2[$i]\t\t$num2[$i]\n");
}

#sub1
########
sub sum
########
{
	my $sum_cur = 0;
	foreach(@_){
		$sum_cur += $_;
	}
	$sum_cur;
}
#sub2
###############
sub array_cout
###############
{
	my $cout = 0;
	my $tar_cout = shift;
	foreach(@_){
		($_ < $tar_cout) ? ($cout++) : return($cout-1);
	}
	return($cout);
}
#sub3
################
sub single_stat
################
{
	#single_stat($inf,\@lv,\@size,\@num,$lenf,$cutl,$becon,$nnum)

	my $inf = $_[0];
	my @lv = @{$_[1]};
    my $cutl = ($_[5] || 100);
    my $becon = ($_[6] || 0);
    !$becon && (`head -1 $inf` =~ /\blength=*\s*\d+/) && ($becon = 1);
	my @leng;
    get_leng($inf,$_[4],$becon,$cutl,\@leng,$_[7]);#sub3.1
	my $toal_leng = sum(@leng);#sub1
	my $toal_num = @leng;
	my $longest = $leng[-1];
	my (@size,@num);
	my ($i,$cur_leng) = (-1, 0);	
	foreach my $rate(@lv){
		my $tar_leng = $toal_leng * (100 - $rate) /100;
		until(($cur_leng > $tar_leng) || ($i == $#leng)){
			$i++;
			$cur_leng += $leng[$i];
		}
		push @size,$leng[$i];
		push @num,($toal_num - $i);
	}
    my $cout_k1 = ($cutl < 100) ? ($toal_num - array_cout(100,@leng)) : $toal_num;
	my $cout_k2 = $toal_num - array_cout(2000,@leng);
	@{$_[2]} = (@size,$longest,$toal_leng,'-','-');
	@{$_[3]} = (@num,'-','-',$cout_k1, $cout_k2);
}
#sub3.1
#############
sub get_leng
#############
{
    my ($infa,$outlen,$becon,$cutl,$leng,$nnum) = @_;
    if($becon){
       my $get_conl;
       if($outlen && !(-s $outlen)){
           $get_conl = 'perl -ne \'/^>(\S+)/ || next;$id=$1;/length=*\s*(\d+)/ || next;'. "(\$1<$cutl) && next;" . 'print "$id\t$1\n";\'';
          `$get_conl $infa > $outlen`;
       }else{
           $get_conl = 'perl -ne \'/^>(\S+)/ || next;$id=$1;/length=*\s*(\d+)/ || next;'. "(\$1<$cutl) && next;" . 'print "$1\n";\'';
           chomp(@{$leng} = `$get_conl $infa`);
           @{$leng} = sort {$a<=>$b} @{$leng};
           return(0);
       }
    }
    if($outlen && (-s $outlen)){
        chomp(@{$leng} = `awk \'(\$2>=$cutl){print \$2}\' $outlen`);
        @{$leng} = sort {$a<=>$b} @{$leng};
        return(0);
    }
    open INFA,$infa;
    $outlen && !(-s $outlen) && (open OUTLEN,">$outlen");
    $/=">";<INFA>;$/="\n";
    my $out;
    while(<INFA>){
        my($id,$seq,$len,$num);
        /^(\S+)/ && ($id = $1);
        $/=">";chomp($seq = <INFA>);$/="\n";
        $seq =~ s/\s+//g;
        $len = length($seq);
        ($len < $cutl) && next;
        $leng && (push @{$leng},$len);
        (!$outlen || (-s $outlen)) || next;
        if($nnum){
            $num = (($seq=~s/N//ig) || 0);
            $out .= "$id\t$len\t$num\n";
        }else{
            $out .= "$id\t$len\n";
        }
    }
    close INFA;
    $outlen && !(-s $outlen) && ((print OUTLEN $out),(close OUTLEN));
    $leng && (@{$leng} = sort {$a<=>$b} @{$leng});
}
