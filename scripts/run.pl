#!/opt/local/bin/perl
# perl script to run maggs-rossetto code
# then do some bookkeeping with the output
#
# Reads in the input file for the fortran program and also the
# command-line arguments provided to it, then generates input files
# for every combination thereof and executes them one by one.
# Most useful for submitting batches of jobs to the HPC cluster,
# hence the hardcoded job script halfway down which is also generated.

use strict;
use warnings;
use Env;
use Cwd;
use Getopt::Long;
use File::Path qw(make_path);
use File::Copy;
use File::Basename;
use Data::Dumper qw(Dumper);

my $row;
my $fh;
my $fh_temp;
my $input;
my $help = '';
my $dorun = 1;
my $doplots = 1;
my $doanalysis = 1;
my $docontour = 1;
my $doquiver = 1;
my $dolorentz = 1;
my $doquadrics = 1;
my $dohelmholtz = 1;
my $arrow_width = 0.002;
my $lgf_path_fix = 0;
my $dpi = 200;
my $inputfile = 'in/start.in';
my $tempinputfile = '';
my $comment = '';
my @lengths = '';
my @temperatures = '';
my @charges = '';
my @charge_values = '';
my @core_energies = '';
my @spacings = '';
my @params_temp = '';
my @stamparray = '';
my $stampdir = '';

# Get command line options; they all have (hopefully) sensible defaults
$input = GetOptions ("help"=> \$help,
                     "run=i"=> \$dorun,
                     "analysis=i"=> \$doanalysis,
                     "plot=i"=> \$doplots,
                     "contour=i"=> \$docontour,
                     "quiver=i"=> \$doquiver,
                     "quadrics=i"=> \$doquadrics,
                     "helmholtz=i"=> \$dohelmholtz,
                     "lorentz=i"=> \$dolorentz,
                     "directory=s"=> \$stampdir,
                     "lengths=s"=> \@lengths,
                     "temperatures=s"=> \@temperatures,
                     "charges=s"=> \@charges,
                     "charge_values=s"=> \@charge_values,
                     "core_energies=s"=> \@core_energies,
                     "spacings=s"=> \@spacings,
                     "input_file=s"=> \$inputfile,
                     "comment=s"=> \$comment);

# Print help if needed. Exit if something weird's been chucked in.
my $helpstring = "Run script for Maggs-Rossetto code. Options:
  --help or -h: Print this help and die.

  --plot or -p = 0/1: Switch plotting of results on/off. Default is 1, set to 0 to turn plotting off.

  --run or -r = 0/1: Switch actual running of code. Default is 1; set to 0 to test the script without running.

  --temperatures or -t = comma-separated list of temperatures. Specifying will cause the script to run once using each temperature given. If a list is given here, the temperature in the input file will be ignored.

  --charge_values  or -ch = comma-separated list of charge values. Specifying will cause the script to run once using each charge value given (i.e. changing the strength of Coulomb interactions, basically). If a list is given here, the temperature in the input file will be ignored.

  --core_energies  or -cor = comma-separated list of charge values. Specifying will cause the script to run once using each charge value given (i.e. changing the strength of Coulomb interactions, basically). If a list is given here, the temperature in the input file will be ignored.

  --spacings or -s = comma-separated list of spacings. Specifying will cause the script to run once using each spacing given. Here for completeness really. If a list is given here, the temperature in the input file will be ignored.

  --input_file or -i = \$filename: Give input file name. Default is \$PWD/in/start.in.

  --comment or -com = \$string: Add a comment to the JSON string I output with each run. I basically use this as a memory aid.
  ";

if ($help || !$input) {
  print $helpstring;
  exit $input;
}

# note: could make an array of references to these and then
# iterate over that, would clean up a couple of bits here
if (@lengths) {
  @lengths = split(/,/,join(',',@lengths));
  shift(@lengths);
}

if (@temperatures) {
  @temperatures = split(/,/,join(',',@temperatures));
  shift(@temperatures);
}

if (@charges) {
  @charges = split(/,/,join(',',@charges));
  shift(@charges);
}

if (@charge_values) {
  @charge_values = split(/,/,join(',',@charge_values));
  shift(@charge_values);
}

if (@core_energies) {
  @core_energies = split(/,/,join(',',@core_energies));
  shift(@core_energies);
}

if (@spacings) {
  @spacings = split(/,/,join(',',@spacings));
  shift(@spacings);
}

# useful directories to keep track of
my $basedir = getcwd();
my $outdir = "$basedir/out";
my $logdir = "$outdir/logs";
print "Creating directory $logdir\n";
make_path($logdir);


my %parameters = get_parameters("$inputfile");

# ensures that none of the lists will be empty
push @lengths, $parameters{L};
push @temperatures, $parameters{temperature};
push @charges, $parameters{charges};
push @charge_values, $parameters{charge_value};
push @core_energies, $parameters{e_c};
push @spacings, $parameters{lattice_spacing};

# ensures we don't waste time doing identical simulations
@lengths = uniq(@lengths);
@temperatures = uniq(@temperatures);
@charges = uniq(@charges);
@charge_values = uniq(@charge_values);
@core_energies = uniq(@core_energies);
@spacings = uniq(@spacings);

# portably change relative path names into absolute ones
my $progname = get_progname("EXEC");

# one line in each timestamped folder with parameters
for( my $i = 0; $i < @temperatures; $i++) {
  for( my $j = 0; $j < @charge_values; $j++) {
    for( my $k = 0; $k < @spacings; $k++) {
      for( my $l = 0; $l < @charges; $l++) {
        for( my $m = 0; $m < @lengths; $m++) {
          for( my $n = 0; $n < @core_energies; $n++) {

            $parameters{temperature} = $temperatures[$i];
            $parameters{charge_value} = $charge_values[$j];
            $parameters{lattice_spacing} = $spacings[$k];
            $parameters{charges} = $charges[$l];
            $parameters{L} = $lengths[$m];
            $parameters{e_c} = $core_energies[$n];
            print "New temperature: $parameters{temperature}\n";
            print "New number of charges: $parameters{charges}\n";
            print "New charge value: $parameters{charge_value}\n";
            print "New core energy: $parameters{e_c}\n";
            print "New lattice spacing $parameters{lattice_spacing}\n";
            print "New system size $parameters{L}\n";

            my $timestamp = get_timestamp();

            # don't think i can call a function inside an array assignment
            # my @stamparray = ($timestamp,'L',$parameters{L},'T', $temperatures[$i],'chg', $charges[$l],'q', $charge_values[$j],'a', $spacings[$k]);
            if ($parameters{canon} =~ /T/ || $parameters{canon} =~ /Y/) {
              @stamparray = ('ce','T', $temperatures[$i],'chg', $charges[$l],$comment);
            } else {
              # gonna want to add core-energy in here probably
              @stamparray = ('gce','T', $temperatures[$i],'e_c',$core_energies[$n],$comment);
            }
            my $stamp = join('_', @stamparray);
            $stampdir = "$outdir/$stamp";
            print "Creating directory $stampdir .\n";
            make_path($stampdir);
            $parameters{stamp} = "$stamp";
            # my $parameters_json = create_json(%parameters);

            if ($doplots) {
              print "Creating directory $stampdir/plots .\n";
              make_path("$stampdir/plots");
            }

            my $logfile = "$logdir/$stamp.log";

            foreach (keys %parameters) {

              if ($lgf_path_fix = 0) {
                if ($_ =~ /lgf_path/) {
                  # my ($fnm,$dirs,$suff) = fileparse($parameters{$_});
                  $parameters{$_} = "$basedir/$parameters{$_}";
                  $lgf_path_fix = 1;
                }
              }

              if ($_ =~ /_file/) {
                my ($fnm,$dirs,$suff) = fileparse($parameters{$_});
                $parameters{$_} = "$stampdir/$fnm$suff";
              }
            }

            $tempinputfile = "$basedir/$stamp.temp";

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

            # generate job file
            my $genjobfile = 0;
            if ($genjobfile) {
              my $email = q(zcapg55@ucl.ac.uk);
              my $vf = '960M';
              my $pe = '32';
              my $jobname = "mr_T_$parameters{temperature}_L_$parameters{L}_canon";
              my $jobfilecontent = qq(
#!/bin/bash -f
# ------------------------------
#\$ \-M $email
#\$ \-m bes
#\$ \-V
#\$ \-j y
#\$ \-cwd
#\$ \-N '$jobname'
#\$ \-S /bin/bash
#\$ \-l vf=$vf
#\$ \-pe ompi $pe
#
echo "Got \${NSLOTS} slots."
IPWD=`pwd`
echo "in \${IPWD}."
mpirun --mca btl ^openib --mca mtl ^psm --n \${NSLOTS} $basedir/$progname -v $tempinputfile
exit 0
);
              my $jobfilename = "$basedir/$stamp.job";
              write_to_file($jobfilename, $jobfilecontent, "write");

              my $runcmd = qq(qsub -q 32-32.q $jobfilename);
              print "Running $progname with command $runcmd\n";
              if ($dorun) {
                system($runcmd);
              }

            } else {

              # run program
              my $np = 1;
              # my $mpi_args = '--bind-to none';
              my $mpi_args = '';
              my $verbose = '-v';
              my $runcmd = qq(mpirun -np $np $mpi_args $progname $verbose $tempinputfile 2>&1 | tee $logfile);
              print "Running $progname with command $runcmd\n";
              if ($dorun) {
                system($runcmd);
              }

            } # end of job file???

            copy($tempinputfile, "$stampdir/input.in");
            # unlink($tempinputfile);
            print "Do analysis = $doanalysis\n";
            if ($doanalysis) {
              print "Doing analysis\n";
              my $analysis_script = "$basedir/scripts/analyse.pl";
              system(qq[$analysis_script -d=$stampdir --plot=$doplots --contour=$docontour --quiver=$doquiver --lorentz=$dolorentz --quadrics=$doquadrics --helmholtz=$dohelmholtz]);
            }

            
          }
        }
      }
    }
  }
}

sub get_timestamp {

  # Copying everything to timestamped directories ensures I keep backups
  (my $sec, my $min, my $hour, my $day , my $monthoffset , my $yearoffset , my $weekday , my $yearday , my $dst) = localtime();
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

  return $timestamp;

}

sub get_progname {

  # searches for the string you give it in the makefile. i use "EXEC"
  my ($progvar) = @_;
  my $progname = qx[awk '/^$progvar*/ { print \$3 }' makefile];
  print "Name of executable: $progname";
  chomp($progname);

  return $progname;

}

sub get_parameters {

  # get parameters from input file; add them to JSON file
  my ($input) = @_;
  open my $fh, '<:encoding(UTF-8)', "$input"
    or die "Unable to open file:$!\n";
  my %parameters = map { split /\s+/; } <$fh>;
  close $fh;

  $parameters{comment} = "$comment";

  if (!@temperatures) {
    push @temperatures, $parameters{temperature};
  }
  if (!@charges) {
    push @charges, $parameters{charges};
  }
  if (!@charge_values) {
    push @charge_values, $parameters{charge_value};
  }
  if (!@spacings) {
    push @spacings, $parameters{lattice_spacing};
  }

  return %parameters;

}

sub write_to_file {

  # writes, not appends!!!
  # print "\nwrite_to_file\n";
  @params_temp = @_;
  if ($params_temp[2] =~ /write/) {
    open($fh_temp, '>:encoding(UTF-8)', $params_temp[0])
    or die "Unable to open file $params_temp[0]\n";
  } elsif ($params_temp[2] =~ /append/) {
    open($fh_temp, '>>:encoding(UTF-8)', $params_temp[0])
    or die "Unable to open file $params_temp[0]\n";
  } else {
    die "Bad argument passed to write_to_file.\n";
  }
  print $fh_temp "$params_temp[1]";
  close $fh_temp;
  return;

}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
