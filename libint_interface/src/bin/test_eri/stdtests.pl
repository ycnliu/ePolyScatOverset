#!/usr/bin/perl

#
# collection of standard ERI tests
# tests most features of the generator
#

use Getopt::Long;

# get optional command-line arguments
my $runall = 0;
my $runlong = 0;
&GetOptions("all" => \$runall,
	    "long" => \$runlong);
$runlong = 1 if $runall == 1;

(usage() and die) if ($#ARGV > 0);

my $tester = "test_eri.pl";
my $timer = "time_eri.pl";

my @tester_args = ("1 0 1 0 0", "1 0 1 0 10000", "1 1 1 1 0", "1 1 1 1 10000", "2 0 2 0 0", "2 1 2 0 10000", "1 1 1 1 10000 16 0", "1 1 1 1 10000 16 1", "1 1 1 1 10000 1 1", "1 1 1 1 0 16 0", "1 1 1 1 0 16 1", "1 1 1 1 10000 16 0 1") ;
my @tester_args_long = ("3 0 3 0 10000", "4 0 4 0 10000", "2 2 2 2 0", "2 2 2 2 10000") ;

if ($runall || !$runlong) {
  foreach my $tester_arg (@tester_args) {
    printf STDOUT "Test \"$tester_arg\" ... ";
    system("./$tester $tester_arg 1>/dev/null 2>/dev/null") && die("failed");
    printf STDOUT "passed\n";
    system("make genclean 1>/dev/null 2>/dev/null");
  }
}

if ($runlong) {
  printf STDOUT "Starting time-consuming tests\n";
  foreach my $tester_arg (@tester_args_long) {
    printf STDOUT "Test \"$tester_arg\" ... ";
    system("./$tester $tester_arg 1>/dev/null 2>/dev/null") && die("failed");
    printf STDOUT "passed\n";
    system("make genclean 1>/dev/null 2>/dev/null");
  }
}

sub usage {
  printf STDERR "USAGE: stdtest.pl [--all] [--long]\n";
  printf STDERR "       Options:\n";
  printf STDERR "         --all              -- run time-consuming tests also (by default, only short tests are attempted)\n";
  printf STDERR "         --long             -- run time-consuming tests only\n";
}
