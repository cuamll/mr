#!/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec;

my $sec;
my $min;
my $hour;
my $day;
my $monthoffset;
my $yearoffset;
my $weekday;
my $yearday;
my $dst;

($sec, $min, $hour, $day, $monthoffset, $yearoffset, $weekday, $yearday, $dst) = localtime();
my $year = 1900 + $yearoffset;
my $month = 1 + $monthoffset;
my $sep = '_';
my $timestamp = join $sep, $year, $month, $day, '', $hour, $min;
print $timestamp . "\n";
my $basedir = File::Spec->curdir();

# stdout will be copied to logs, with timestamp
# newest set of files get copied to directory with matching timestamp
my $logdir = "$basedir/out/logs";
my $timedir = "$basedir/out/$timestamp";
make_path($logdir,$timedir);

my $logfile = $logdir . "/" . $timestamp . ".log";

my $progname = "mr_test";
my $inputfile = "$basedir/in/start.in";

my $runcmd = qq(./$progname $inputfile 2>&1 | tee $logfile);
print $runcmd;

#system($runcmd);
