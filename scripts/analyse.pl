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
my $domultiplot = 1;
my $dorun = 1;
my $dospr = 1;
my $doquiver = 1;
my $dolorentz = 1;
my $doquadrics = 1;
my $dohelmholtz = 1;
my $arrow_width = 0.002;
my $core_energy = 0.0;
my $dpi = 200;
my $inputfile = '';
my $comment = '';
my $stampdir = '';

# Get command line options; they all have (hopefully) sensible defaults
$input = GetOptions ("help"=> \$help,
                     "plot=i"=> \$doplots,
                     "spr=i"=> \$dospr,
                     "quiver=i"=> \$doquiver,
                     "multiplot=i"=> \$domultiplot,
                     "lorentz=i"=> \$dolorentz,
                     "quadrics=i"=> \$doquadrics,
                     "helmholtz=i"=> \$dohelmholtz,
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

# the lines corresponding to the fourier space helmholtz
# decomposition probably aren't there: add them
my %extra_params = ();
$extra_params{"s_ab_t_file"} = "$stampdir/s_ab_t.dat";
$extra_params{"s_ab_l_file"} = "$stampdir/s_ab_l.dat";
$extra_params{"s_perp_t_file"} = "$stampdir/s_perp_t.dat";
$extra_params{"s_perp_l_file"} = "$stampdir/s_perp_l.dat";
$extra_params{"s_par_t_file"} = "$stampdir/s_par_t.dat";
$extra_params{"s_par_l_file"} = "$stampdir/s_par_l.dat";
$extra_params{"charge_generation"} = "DIPOLE";

foreach (keys %extra_params) {
  $parameters{$_} = $extra_params{$_} unless (exists($parameters{$_}));
}

foreach (keys %parameters) {
  if ($_ =~ /_file/) {
    my ($fnm,$dirs,$suff) = fileparse($parameters{$_});
    $parameters{$_} = "$stampdir/$fnm$suff";
  }
}

my $tempinputfile = "$stampdir/temp";

open($fh, '>:encoding(UTF-8)', $tempinputfile)
or die "Unable to create temporary input file:$!\n";
print "Creating temporary input file at $tempinputfile\n";

for my $key (keys %parameters) {
if ($key !~ /comment/) {
  if ($key !~ /stamp/) {
    print $fh "$key $parameters{$key}\n";
  }
}
}
close $fh;

# some of my old results don't include this parameter
if (exists $parameters{'e_c'}) {
  $core_energy = $parameters{'e_c'}
}

# job file tells us how many MPI threads were used in simulation
my ($filename,$directories,$suffix) = fileparse($inputfile);
my @arr = split(/\//,$directories);
my $jobfile = "$directories$arr[-1].job";
print "Job file = $jobfile\n";

# grab the number of MPI slots used. Bit of a hack but ¯\_(ツ)_/¯
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
  $doquadrics,
  $dospr,
  $doplots,
  $doquiver,
  $domultiplot,
  $dolorentz,
  $dohelmholtz,
);

my @file = (
 "$basedir/scripts/quadrics.py",
 "$basedir/spr_fh",
 "$basedir/scripts/plot.pl",
 "$basedir/scripts/quiver.py",
 "$basedir/scripts/sab_multiplot.py",
 "$basedir/scripts/lorentz.py",
 "$basedir/scripts/helmholtz.py"
);

# NB: for the plotfile especially there are extra parameters;
# I use the defaults here but check the plot file if they need changing
my @cmd = (
  qq[python $file[0] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi],
  qq[$file[1] $tempinputfile],
  qq[$file[2] -d=$stampdir -s="$no_slots"],
  qq[python $file[3] $stampdir $parameters{L} $arrow_width $dpi],
  qq[python $file[4] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi],
  qq[python $file[5] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi],
  qq[python $file[6] $stampdir $parameters{L} $parameters{temperature} $core_energy $dpi],
);

for my $i (0..$#run) {
  if ($run[$i]) {
    print "Running $file[$i]\n";
    system($cmd[$i]);
  }
}

# unlink $tempinputfile;

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
