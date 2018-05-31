#!/opt/local/bin/perl
# perl script to analyse results from maggs-rossetto CG code
use strict;
use warnings;
use Env;
use Cwd;
use Getopt::Long;
use File::Basename;

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
my $core_energy = 0.0;
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
                     "quadrics=i"=> \$doquadrics,
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
print "Input file: $inputfile\n";
my %parameters = get_parameters("$inputfile");

foreach (keys %parameters) {
  if ($_ =~ /_file/) {
    my ($fnm,$dirs,$suff) = fileparse($parameters{$_});
    $parameters{$_} = "$stampdir/$fnm$suff";
  }
}

# some of my old results don't include this parameter
if (exists $parameters{'e_c'}) {
  $core_energy = $parameters{'e_c'}
}

# job file tells us how many MPI threads were used in simulation
my ($filename,$directories,$suffix) = fileparse($inputfile);
my @arr = split(/\//,$directories);
my $jobfile = "$directories/$arr[-1].job";
print "Job file = $jobfile\n";

# grab the number of MPI slots used. Hacky af hehe
if (-f $jobfile) {
  open $fh, '<', $jobfile;
  while ($line = <$fh>) {
    # print $line;
    if ($line =~ /ompi/) {
      $no_slots = substr $line, -3;
      print "Line = $line\nNumber of slots = $no_slots\n";
    }
  }
  close $fh;
}

my @run = (
  $doplots,
  $docontour,
  $doquiver,
  $dolorentz,
  $doquadrics
);

my @file = (
 "$basedir/scripts/plot.pl",
 "$basedir/scripts/s_perp_contours.py",
 "$basedir/scripts/quiver.py",
 "$basedir/scripts/lorentz.py",
 "$basedir/scripts/quadrics.py"
);

# NB: for the plotfile especially there are extra parameters;
# I use the defaults here but check the plot file if they need changing
my @cmd = (
  qq[$file[0] -d=$stampdir -s="$no_slots"],
  qq[python $file[1] $stampdir $parameters{L} $dpi],
  qq[python $file[2] $stampdir $parameters{L} $arrow_width $dpi],
  qq[python $file[3] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi],
  qq[python $file[4] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi]
);

for my $i (0..$#run) {
  if ($run[$i]) {
    print "Running $file[$i]\n";
    system($cmd[$i]);
  }
}

sub get_parameters {

  # get parameters from input file; add them to JSON file
  my ($tempinput) = @_;
  open my $tempfh, '<:encoding(UTF-8)', "$tempinput"
    or die "Unable to open file:$!\n";
  my %tempparameters = map { split /\s+/; } <$tempfh>;
  close $tempfh;

  $tempparameters{comment} = "$comment";

  return %tempparameters;

}
