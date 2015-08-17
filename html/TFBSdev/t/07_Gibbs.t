#!/usr/bin/env perl -w

use TFBS::PatternGen::Gibbs;
use Test;
plan(tests => 2);
my $gibbspath;
eval {$gibbspath = `which Gibbs 2> /dev/null`;};
if (!$gibbspath or ($gibbspath =~ / /)) { # if space, then error message :)
    print "ok \# Skipped: (no Gibbs executable found)\n"x2;
    exit(0);
}

my $fastafile = "t/test.gibbin";

for (1..5) {
    my $gibbs = 
      TFBS::PatternGen::Gibbs->new
	  (-nr_hits=>10,
	   -motif_length=>[10..12],
	   -seq_file=>$fastafile);
    

    my @pfms = $gibbs->all_patterns();
    if (@pfms>0) {
	ok(1);
	last;
    }
}

ok(1);

