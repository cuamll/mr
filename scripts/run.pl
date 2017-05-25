#!/bin/perl
# perl script to run maggs-rossetto code
# then do some bookkeeping with the output
use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec;
use File::Copy;
use File::Basename;
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

# predictable time/date formatting
my @timedate = ($year, $month, $day, $hour, $min);
foreach (@timedate) {
  if ($_ < 10) {
    $_ = "0" . $_;
  }
}
my $sep = '_';
my $timestamp = join($sep, @timedate);
#print "$timestamp\n";

# stdout will be copied to logs, with timestamp
# newest set of files get copied to directory with matching timestamp
my $basedir = File::Spec->curdir();
my $outdir = "$basedir/out";
my $logdir = "$outdir/logs";
my $timedir = "$outdir/$timestamp";
print "Creating directories $logdir and $timedir .\n";
#make_path($logdir,$timedir);

my $logfile = "$logdir/$timestamp.log";

# NB: this only works if your -o variable in makefile starts with EXEC
my $progvar = "EXEC";
my $progname = qx[awk '/^$progvar*/ { print \$3 }' makefile];
print "Name of executable: $progname";
chomp($progname);

# this could be given as an argument, I guess?
my $inputfile = "$basedir/in/start.in";

# get parameters from input file; add them to JSON file
open my $fh, '<:encoding(UTF-8)', "$inputfile"
  or die "Unable to open file:$!\n";
my %parameters = map { split /\s+/; } <$fh>;
close $fh;
$parameters{timestamp} = "$timestamp";
my $parameters_json = encode_json [%parameters];

# ----------
# uncomment this block when done testing

# one line in each timestamped folder with parameters
#open $fh, '>:encoding(UTF-8)', "$timedir/parameters.json"
#  or die "Unable to open JSON file:$!\n";
#print $fh "Adding JSON parameter data to " .
#          "$timedir/parameters.json:\n$parameters_json\n";
#close $fh;

# keep a list of every run and its parameters
#my $jsondb = "$logdir/db.json";
#open $fh, '>>:encoding(UTF-8)', "$jsondb"
#  or die "Unable to open JSON file:$!\n";
## this might not work idk the syntax
#print $fh "$parameters_json\n";
#close $fh;
# ----------

# run program
my $runcmd = qq(./$progname $inputfile 2>&1 | tee $logfile);
print "Running $progname with command $runcmd\n";

#system($runcmd);

# program should have terminated now
# we want to copy the out directory into the timestamped one.
# do ls, make a list, then iterate over the list and copy

opendir(my $dh, $outdir) or die "Can't open output directory: $!";
foreach (keys %parameters) {
  if ($_ =~ /_file/) {
    my $value = %parameters{$_};
    print "$parameters{$_}\n";
    my $file;
    my $dir;
    my $ext;
    ($file, $dir, $ext) = fileparse($value);
    #print $dir . $file .$ext . "\n";
    #copy($dir . $file . $ext, "$timedir/$file" . $ext) or die "Copy failed: $!";
    print "Copying $dir$file$ext to $timedir/$file$ext\n";
    #copy($dir . $file . $ext, "$timedir/$file" . $ext) or die "Copy failed: $!";
  }
}
close($dh);

# then probably averaging stuff on the raw output, e_fields and charge dist.
# then gnuplot
