#!/usr/bin/perl

# yangm@idsse.ac.cn 2023/04/01 01:54
# Updated to use thread pool model with Thread::Queue

use strict;
use warnings;
use threads;
use Thread::Queue;

unless(@ARGV > 0){
    print STDERR "Usage: $0 <wkdir> <postfix_sh> <parallel_number>\n";
    exit 0;
}

my($wkdir, $postfix, $parallel_number) = @ARGV;
$wkdir ||= ".";
$postfix ||= ".sh";
$parallel_number ||= 2;

die "The directory \"$wkdir\" does not exist" unless (-d $wkdir);

# Search for all scripts files
my @scripts = glob("$wkdir/*$postfix");
die "No script files fitted for \"$wkdir/*$postfix\"" unless(@scripts);

# Create task queue
my $task_queue = Thread::Queue->new(@scripts);

# Launch a fixed number of threads, each continuously fetching and executing tasks from a queue.
my @threads;
foreach(1..$parallel_number){
    push @threads, threads->create(
        sub {
            while(defined(my $script = $task_queue->dequeue_nb)){
                print "Thread ", threads->self->tid, " executing $script\n";
                system("sh $script >$script.log 2>$script.error");
            }
        }
    );
}

# Wait for all threads to complete
$_->join for @threads;

print "All task script file with a surfix of \"$postfix\" finished!!!\n";

__END__
