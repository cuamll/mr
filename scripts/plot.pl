#!/opt/local/bin/perl
# PRELIMINARY!
# i'm probably gonna expand the scripting to automate more stuff.
# currently this just plots heatmaps of correlation functions
# call with "perl plot.pl -l=L ..." etc. from scripts dir
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

my $three_d = 0;
my $dir;
my $kz;
my $palette = 'inferno.pal';
my @inputfiles;
my @outputfiles;
my @tempfiles;
my @gnuplotargs;
my $input;
$input = GetOptions ("d=s"=> \$dir,
                     "p=s"=> \$palette);

if ($palette !~ /.*\.pal/) {
  $palette = "$palette.pal";
}

# we pass the absolute path to the timestamp directory. now to get parameters:
my %parameters = get_parameters("$dir/input.in");

my @filenames;
foreach my $key (keys %parameters) {
  if ($key =~ /_file/) {
    if ($parameters{$key} =~ /.*s_.*/) {
      #if ($parameters{$key} !~ /.*s_ab.*/) {
        (my $tempfile, my $tempdir, my $tempext) = fileparse($parameters{$key}, qr/\.[^.]*/);
        push @filenames, $tempfile
      #}
    }
  }
}
my @titles;
my $s_component; my $field_component; my $s_string; my $field_string;
my $linetitle;

foreach my $file (@filenames) {

  if ($file =~ /s_([a-z]+)_([a-z]+)/) {
    $s_component = $1;
    $field_component = $2;

    if ($s_component =~ /xx/) {
      $s_string = qq(S^{$s_component} - );
    } elsif ($s_component =~ /perp/) {
      $s_string = qq(S^{⟂} - );
    } elsif ($s_component =~ /par/) {
      $s_string = qq(S^{∥} - );
    } elsif ($s_component =~ /ab/) {
      # print "S^{/alpha /beta still there???}";
      $s_string = qq(S^{/Symbol a /Symbol b} - );
    } else {
      die "s_component is weird: $s_component $!\n";
    }

    if ($field_component =~ /total/) {
      $field_string = qq(total);
    } elsif ($field_component =~ /irrot/) {
      $field_string = qq(irrotational);
    } elsif ($field_component =~ /rot/) {
      $field_string = qq(rotational);
    } else {
      die "field_component is wrong: $field_component $!\n";
    }

    $linetitle = $s_string . $field_string;
  } elsif ($file =~ /s_([a-z]+)/) {

    if ($1 =~ /charge/) {
      $linetitle = qq(g^{+-}(k));
    } elsif ($1 =~ /direct/) {
      $linetitle = qq(g^{+-}(r));
    } else {
      die "File name doesn't match regex. $file $!\n";
    }

  } else {
    die "File name doesn't match anything. $file $!\n";
  }

  push @titles, $linetitle;
}

my $meas = $parameters{measurement_sweeps} / $parameters{sample_interval};

my $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, charge value = $parameters{charge_value}.\n\n$meas measurements from $parameters{measurement_sweeps} MC steps.);

my $basedir = File::Spec->curdir();
my $plotpath = "$basedir/scripts/";
my $plotsuffix = ".png";
my $gnuplotscript = "$plotpath" . "heatmap.p";

my $inpath = "$dir/";
my $insuffix = ".dat";
my $outpath = "$dir/plots/";
my $tempsuffix = '.temp';

for my $i (0..$#filenames) {
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @tempfiles, $inpath . $filenames[$i] . $tempsuffix;
  push @outputfiles, $outpath . $filenames[$i] . $plotsuffix;

  # don't need any of this in 2d
  # in 3d, create temp files for each slice in kz
  if ($three_d) {

    my $length = $parameters{L};
    my $awkcall;
    my $bz = 2;
    my $sperp_size = 5;
    my $lowerbound = ((($bz * $length) + 1)**2) * ($kz + ($bz * ($length / 2)));
    my $upperbound = $lowerbound + 1 + (($bz * $length) + 1)**2;
    my $splb = ((($sperp_size * $length) + 1)**2) * ($kz + ($sperp_size * ($length / 2)));
    my $spub = $splb + 1 + (($sperp_size * $length) + 1)**2;

    if (index($inputfiles[$i],'s_perp') != -1) {
      $awkcall = qq[awk '{ if (NR>$splb && NR<$spub) print \$1,\$2,\$4 }' $inputfiles[$i] > $tempfiles[$i]];
    } else {
      $awkcall = qq[awk '{ if (NR>$lowerbound && NR<$upperbound) print \$1,\$2,\$4 }' $inputfiles[$i] > $tempfiles[$i]];
    }
    system($awkcall);

  }
}

warn "Different number of titles and files!\n" unless @titles == @filenames;

for my $i (0..$#filenames) {
  my $syscall;

  if ($three_d) {
    push @gnuplotargs, qq(FILE='$tempfiles[$i]'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$plottitle'; LINETITLE = '$titles[$i]'; PALETTE = '$palette');
  } else {

    push @gnuplotargs, qq(FILE='$inputfiles[$i]'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$plottitle'; LINETITLE = '$titles[$i]'; PALETTE = '$palette';);

    #if (index($inputfiles[$i],'s_direct') != -1) {
    #  push @gnuplotargs, qq( PITICS = 'N');
    #} else {
    #  push @gnuplotargs, qq( PITICS = 'Y');
    #}

  }

  $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  system($syscall);

  # delete the temp files
  unlink($tempfiles[$i]);
}

sub get_parameters {

  # get parameters from input file; add them to JSON file
  my ($input) = @_;
  open my $fh, '<:encoding(UTF-8)', "$input"
    or die "Unable to open file:$!\n";
  my %parameters = map { split /\s+/; } <$fh>;
  close $fh;

  return %parameters;

}
