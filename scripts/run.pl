#!/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec;
use JSON::MaybeXS qw(encode_json decode_json);
use Data::Dumper qw(Dumper);

my $sec;
my $min;
my $hour;
my $day;
my $monthoffset;
my $yearoffset;
my $weekday;
my $yearday;
my $dst;
my $row;

($sec, $min, $hour, $day, $monthoffset, $yearoffset, $weekday, $yearday, $dst) = localtime();
my $year = 1900 + $yearoffset;
my $month = 1 + $monthoffset;

my @timedate = ($year, $month, $day, $hour, $min);
foreach (@timedate) {
  if ($_ < 10) {
    $_ = "0" . $_;
  }
}
my $sep = '_';
my $timestamp = join($sep, @timedate);
print "$timestamp\n";
my $basedir = File::Spec->curdir();

# stdout will be copied to logs, with timestamp
# newest set of files get copied to directory with matching timestamp
my $logdir = "$basedir/out/logs";
my $timedir = "$basedir/out/$timestamp";
make_path($logdir,$timedir);

my $logfile = $logdir . "/" . $timestamp . ".log";

my $progname = "mr_test";
my $inputfile = "$basedir/in/start.in";

# get parameters from input file; add them to JSON file
open my $fh, '<:encoding(UTF-8)', "$inputfile"
  or die "Unable to open file:$!\n";
my %parameters = map { split /\s+/; } <$fh>;
close $fh;

$parameters{timestamp} = "$timestamp";

#foreach my $name (sort { lc $a cmp lc $b } keys %parameters) {
#  print "$name, $parameters{$name}\n";
#}

print Dumper \%parameters;

my $parameters_json = encode_json [%parameters];

#if (open($fh, '<:encoding(UTF-8)', $inputfile)) {
#  while ($row = <$fh>) {
#    my $key;
#    my $value;
#
#    chomp $row;
#    print $row;
#    my @temp = split /\s+/, $row;
#    #push %parameters, split /\s+/, $row;
#    
#  } 
#} else {
#  warn "Could not open $inputfile. Check path and try again. $!"
#}

my $runcmd = qq(./$progname $inputfile 2>&1 | tee $logfile);
print $runcmd;

#system($runcmd);


