#!/opt/local/bin/perl
# perl script to run maggs-rossetto code
# then do some bookkeeping with the output
use strict;
use warnings;
use Env;
use Cwd;
use Getopt::Long;
use File::Path qw(make_path);
use File::Copy;
use File::Basename;
# use JSON::MaybeXS qw(encode_json decode_json);
use Data::Dumper qw(Dumper);

my $row;
my $fh;
my $line;
my $no_slots;
my $fh_temp;
my $input;
my $help = '';
my $doplots = 1;
my $dorun = 1;
my $docontour = 1;
my $doquiver = 1;
my $dolorentz = 1;
my $doquadrics = 1;
my $arrow_width = 0.002;
my $dpi = 200;
my $inputfile = '';
my $comment = '';
my $stampdir = '';

# Get command line options; they all have (hopefully) sensible defaults
$input = GetOptions ("help"=> \$help,
                     "plot=i"=> \$doplots,
                     "contour=i"=> \$docontour,
                     "quiver=i"=> \$doquiver,
                     "lorentz=i"=> \$dolorentz,
                     "directory=s"=> \$stampdir,
                     "input_file=s"=> \$inputfile,
                     "comment=s"=> \$comment);

# Print help if needed. Exit if something weird's been chucked in.
my $helpstring = "Run script for Maggs-Rossetto code. Options:
  --help or -h: Print this help and die.

  --plot or -p = 0/1: Switch plotting of results on/off. Default is 1, set to 0 to turn plotting off.

  --directory or -d: Give the directory where the relevant output files are.

  --input_file or -i = \$filename: Give input file name. Default is \$PWD/in/start.in.

  --comment or -com = \$string: Add a comment to the JSON string I output with each run. I basically use this as a memory aid.
  ";

if ($help || !$input) {
  print $helpstring;
  exit $input;
}

my $basedir = getcwd();
$inputfile = "$stampdir/input.in";
print $inputfile;
my %parameters = get_parameters("$inputfile");

# get number of slots
my ($filename,$directories,$suffix) = fileparse($inputfile);
my @arr = split(/\//,$directories);
my $jobfile = "$directories/$arr[-1].job";
print "job file = $jobfile\n";

if (-f $jobfile) {
  open $fh, '<', $jobfile;
  while ($line = <$fh>) {
    if ($line =~ /ompi/) {
      $no_slots = substr $line, -3;
      print "Line = $line\nNumber of slots = $no_slots\n";
    }
  }
}

foreach (keys %parameters) {
  if ($_ =~ /_file/) {
    my ($fnm,$dirs,$suff) = fileparse($parameters{$_});
    $parameters{$_} = "$stampdir/$fnm$suff";
  }
}

if ($doplots) {
  print "Creating directory $stampdir/plots .\n";
  make_path("$stampdir/plots");
  my $plotfile = "$basedir/scripts/plot.pl";
  my $measurements = $parameters{measurement_sweeps} / $parameters{sample_interval};
  my $kz = 0;
  my $palette = "~/.config/gnuplot/inferno.pal";
  my $plotcmd = qq[$plotfile -d=$stampdir -p="$palette" -s="$no_slots"];
  system($plotcmd);
}

if ($docontour) {
  my $contourfile = "$basedir/scripts/s_perp_contours.py";
  my $contourcmd = qq[python $contourfile $stampdir $parameters{L} $dpi];
  system($contourcmd);
}

if ($doquiver) {
  my $quiverfile = "$basedir/scripts/quiver.py";
  my $quivercmd = qq[python $quiverfile $stampdir $parameters{L} $arrow_width $dpi];
  system($quivercmd);
}

if ($dolorentz) {
  my $lorentzfile = "$basedir/scripts/fits.py";
  my $lorentzcmd = qq[python $lorentzfile $stampdir $parameters{L} $parameters{temperature} $parameters{e_c} $dpi];
  system($lorentzcmd);
}

if ($doquadrics) {
  my $quadricsfile = "$basedir/scripts/quadrics.py";
  my $quadricscmd = qq[python $quadricsfile $stampdir $parameters{L} $parameters{temperature} $parameters{e_c} $dpi];
  system($quadricscmd);
}

sub get_parameters {

  # get parameters from input file; add them to JSON file
  my ($input) = @_;
  open my $fh, '<:encoding(UTF-8)', "$input"
    or die "Unable to open file:$!\n";
  my %parameters = map { split /\s+/; } <$fh>;
  close $fh;

  $parameters{comment} = "$comment";

  return %parameters;

}
