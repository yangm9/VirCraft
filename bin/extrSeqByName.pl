#!/usr/bin/env perl
#根据序列名称提取FastA文件中的序列并输出FASTA文件
#用法：perl extrSeqByName.pl <序列名称列表> <输入FA文件> <结果FA文件>
use strict;

die "$0 <seq_name_list> <fasta_file> <extracted_fasta> [0/1]\n" if(@ARGV<3);

$ARGV[3]||=0;
#获取序列名称
open LIST,$ARGV[0] or die $!;
my %SeqName=();
while(<LIST>){
    chomp;
    $SeqName{$_}=1;
}
close LIST;

#根据序列名提取FastA序列
$/='>';
open IN,$ARGV[1] or die $!;
open OUT,">$ARGV[2]" or die $!;
<IN>;
my $n=0;
while(<IN>){
    chomp;
    my($id,$seq)=split /\n/,$_,2;
    $seq=~s/\s+//g;
    my $mother_id=(split /\s/,$id,2)[0];
    $mother_id=~s/_\d+$// if($ARGV[3]); #如果添加第四个非0参数，则本行用于去除prodigal的编号
    if(exists $SeqName{$mother_id}){
        print OUT ">$id\n$seq\n";
        #delete $SeqName{$id};
    }
}
$/='\n';
close IN;
close OUT;

#输出不存在序列的序列名称列表
#print "The sequence name without sequence are as follows:";
#foreach my $key(keys %SeqName){
#    print "$key\n";
#}
__END__
