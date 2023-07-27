#!/usr/bin/env perl
#yangm@idsse.ac.cn
use strict;
$/=">";
open FA,$ARGV[0] or die $!;
my %fa=();
<FA>;
my $n=0;
chomp($_=<FA>);
my($info,$seq)=split /\n/,$_,2;
$seq=~s/\s+//g;
my $seqlen=length($seq);
my $name=(split /\s+/,$info,2)[0];
$fa{$name}=$seq;
while(<FA>){
    chomp;
    my($info,$seq)=split /\n/,$_,2;
    my $name=(split /\s+/,$info,2)[0];
    $seq=~s/\s+//g;
    my $nextlen=length($seq);
    die "Some sequence have different lengths!" unless($nextlen==$seqlen);
    $fa{$name}=$seq;
}
$/="\n";
close FA;

my $n=keys %fa;
open PHY,">$ARGV[1]" or die $!;
print PHY "$n $seqlen\n";
foreach my $name(keys %fa){
    my $spaces=&mkSpaces($name);
    print PHY $name.$spaces.$fa{$name}."\n";
}
close PHY;

sub mkSpaces{
    my $name=shift;
    my $num=25-length($name);
    my $spaces=" "x$num;
    return $spaces;
}
