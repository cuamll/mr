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
use JSON::MaybeXS qw(encode_json decode_json);
use Data::Dumper qw(Dumper);

my $row;
my $fh;
my $fh_temp;
my $input;
my $help = '';
my $doplots = 1;
my $dorun = 1;
my $inputfile = 'in/start.in';
my $tempinputfile = '';
my $comment = '';
my @temperatures = '';
my @charges = '';
my @charge_values = '';
my @spacings = '';
my @params_temp = '';

# Get command line options; they all have (hopefully) sensible defaults
$input = GetOptions ("help"=> \$help,
                     "plot=i"=> \$doplots,
                     "run=i"=> \$dorun,
                     "temperatures=s"=> \@temperatures,
                     "charges=s"=> \@charges,
                     "charge_values=s"=> \@charge_values,
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
  
  --spacings or -s = comma-separated list of spacings. Specifying will cause the script to run once using each spacing given. Here for completeness really. If a list is given here, the temperature in the input file will be ignored.

  --input_file or -i = \$filename: Give input file name. Default is \$PWD/in/start.in.

  --comment or -co = \$string: Add a comment to the JSON string I output with each run. I basically use this as a memory aid.
  ";

if ($help || !$input) {
  print $helpstring;
  exit $input;
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

if (@spacings) {
  @spacings = split(/,/,join(',',@spacings));
  shift(@spacings);
}

#foreach (@temperatures) {
#  print "$_\n";
#}


# useful directories to keep track of
my $basedir = getcwd();
my $outdir = "$basedir/out";
my $logdir = "$outdir/logs";
print "Creating directory $logdir\n";
make_path($logdir);

# keep a timestamped copy of output files
$tempinputfile = "$basedir/in.temp";

my %parameters = get_parameters("$inputfile");
print Dumper %parameters;

# portably change relative path names into absolute ones
foreach (keys %parameters) {
  if ($_ =~ /_file/) {
    my $value = $parameters{$_};
    $parameters{$_} = "$basedir/$value";
  }
}

my $progname = get_progname("EXEC");

# one line in each timestamped folder with parameters
if ($dorun) {
  for( my $i = 0; $i < @temperatures; $i++) {
    for( my $j = 0; $j < @charge_values; $j++) {
      for( my $k = 0; $k < @spacings; $k++) {
        for( my $l = 0; $l < @charges; $l++) {

          $parameters{temperature} = $temperatures[$i];
          $parameters{charge_value} = $charge_values[$j];
          $parameters{lattice_spacing} = $spacings[$k];
          $parameters{charges} = $charges[$l];
          print "New temperature: $parameters{temperature}\n";
          print "New number of charges: $parameters{charges}\n";
          print "New charge value: $parameters{charge_value}\n";
          print "New lattice spacing $parameters{lattice_spacing}\n";

          my $timestamp = get_timestamp();

          # don't think i can call a function inside an array assignment
          my @stamparray = ('T', $temperatures[$i],'chg', $charges[$l],'z', $charge_values[$j],'a', $spacings[$k], $timestamp);
          my $stamp = join('_', @stamparray);
          my $stampdir = "$outdir/$stamp";
          print "Creating directories $stampdir and $stampdir/plots .\n";
          make_path($stampdir,"$stampdir/plots");
          $parameters{stamp} = "$stamp";
          my $parameters_json = create_json(%parameters);

          my $logfile = "$logdir/$timestamp.log";

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

          write_to_file("$stampdir/parameters.json", "$parameters_json\n", "write");
          # keep a list of every run and its parameters
          my $jsondb = "$logdir/db.json";
          write_to_file($jsondb, "$parameters_json\n", "append");

          # run program
          my $runcmd = qq(./$progname $tempinputfile 2>&1 | tee $logfile);
          print "Running $progname with command $runcmd\n";
          system($runcmd);

          copy($tempinputfile, "$stampdir/input.in");
          unlink($tempinputfile);

          # program should have terminated now
          # we want to copy the out directory into the timestamped one.
          # do ls, make a list, then iterate over the list and copy
          opendir(my $dh, $outdir) or die "Can't open output directory: $!";
          foreach (keys %parameters) {
            if ($_ =~ /_file/) {
              my $value = $parameters{$_};
              (my $file, my $dir, my $ext) = fileparse($value);
              print "Copying $dir$file$ext to $stampdir/$file$ext\n";
              copy($dir . $file . $ext, "$stampdir/$file" . $ext) or die "Copy failed: $!";
            }
          }
          close($dh);

          # then probably averaging stuff on the raw output, e_fields and charge dist.
          # then gnuplot
          if ($doplots) {
            my $plotfile = "$basedir/scripts/plot.pl";
            my $measurements = $parameters{measurement_sweeps} / $parameters{sample_interval};
            my $kz = 0;
            my $palette = "inferno.pal";
            #my $plotcmd = qq[$plotfile -l=$parameters{L} -t=$parameters{temperature} -m=$measurements -s=$parameters{measurement_sweeps} -c=$parameters{charges} -k=$kz -fp="$stamp" -o="$stampdir/plots/" -p="$palette"];
            my $plotcmd = qq[$plotfile -d=$stampdir -p="$palette"];
            system($plotcmd);
          }
        }
      }
    }
  }
}

########

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

sub create_json {

  # print "\ncreate_json\n";
  my %input = @_;
  my $json = encode_json[%input];
  return $json;

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
