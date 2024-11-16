#!/usr/bin/env perl
# Extracts the sequence from the FastA file based on the sequence name and outputs the FASTA file
# Usage：perl extrSeqByName.pl <sequence_name_list> <input_fasta> <output_fasta>
# yangm@idsse.ac.cn 2022-08-12 18:56 
use strict;

die "$0 <seq_name_list> <fasta_file> <extracted_fasta> [0/1]\n" if(@ARGV<3);

$ARGV[3]||=0;
# Get all sequence name from sequence name list file
open LIST,$ARGV[0] or die $!;
my %SeqName=();
while(<LIST>){
    chomp;
    $SeqName{$_}=1;
}
close LIST;

# Extracts the sequence by sequence name list
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
    $mother_id=~s/_\d+$// if($ARGV[3]); # If a fourth non-0 parameter is added, this line is used to remove prodigal's number
    if(exists $SeqName{$mother_id}){
        print OUT ">$id\n$seq\n";
        delete $SeqName{$id};
    }
}
foreach my $id(keys %SeqName){
    print "$id\n";
}
$/='\n';
close IN;
close OUT;

#Outputs a list of sequence names for which no sequence exists
#print "The sequence name without sequence are as follows:";
#foreach my $key(keys %SeqName){
#    print "$key\n";
#}
__END__
