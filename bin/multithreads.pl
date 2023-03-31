#!/usr/bin/perl
#yangm@idsse.ac.cn 2023/04/01 01:54
use strict;
use threads;
use threads::shared;

unless(@ARGV>0){
    print STDERR "Usage: $0 <wkdir> <postfix_sh> <parallel_number>\n";
    exit 0;
}

my ($wkdir,$postfix,$parallel_number)=@ARGV;
$wkdir||=".";
$postfix||=".sh";
$parallel_number=2;

my @scripts=glob("$wkdir/*$postfix");
my $subsets_ref=&get_subset(\@scripts,$parallel_number);

for my $subset(@$subsets_ref){
    &run_scripts($subset);
}
print "All Done!!!\n";

sub run_scripts($){
    my ($scripts_ref)=@_;
    my @scripts=@$scripts_ref;
    my @threads;
    
    foreach my $script(@scripts) {
        push(@threads, threads->create(sub {
            system("sh $script >$script.o 2>$script.e");
        }));
    }
    foreach my $thread (@threads) {
        $thread->join();
    }
    return 0;
}

sub get_subset($$){
    my ($array_ref, $n) = @_;
    my @array = @$array_ref;
    my @result;
    while (@array) {
        my @subset = splice(@array, 0, $n);
        push(@result, \@subset);
    }
    return \@result;
}
__END__
