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
use Data::Dumper qw(Dumper);

my $three_d = 0;
my $dir; my $kz;
my $palette = '~/.config/gnuplot/jet.pal';
my $addtitles = 1;
my @inputfiles; my @outputfiles; my @tempfiles;
my @gnuplotargs; my @latexargs; my @dvipsargs; my @ps2pdfargs;
my $input = GetOptions ("d=s"=> \$dir,
                     "t=i"=> \$addtitles,
                     "p=s"=> \$palette);

if ($palette !~ /.*\.pal/) {
  $palette = "$palette.pal";
}

# we pass the absolute path to the timestamp directory. now to get parameters:
my %parameters = get_parameters("$dir/input.in");

my $steps = $parameters{no_samples} * $parameters{measurement_sweeps};
my $meas = $steps / $parameters{sample_interval};
my $steps_c = commify($steps);
my $meas_c = commify($meas);

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
my $linetitle; my $plottitle;

if ($addtitles) {
  foreach my $file (@filenames) {

    if ($file =~ /s_([a-z]+)_([a-z]+)/) {
      $s_component = $1;
      $field_component = $2;

      if ($s_component =~ /xx/) {
        $s_string = q($ S^{xx});
      } elsif ($s_component =~ /perp/) {
        $s_string = q($ S^{\perp});
      } elsif ($s_component =~ /par/) {
        $s_string = q($ S^{\parallel});
      } elsif ($s_component =~ /ab/) {
        # print "S^{/alpha /beta still there???}";
        $s_string = q($ S^{\alpha \beta});
      } else {
        die "s_component is wrong: $s_component $!\n";
      }

      if ($field_component =~ /total/) {
        $field_string = q(_{total}$);
      } elsif ($field_component =~ /irrot/) {
        $field_string = q(_{irrotational}$);
      } elsif ($field_component =~ /rot/) {
        $field_string = q(_{rotational}$);
      } else {
        die "field_component is wrong: $field_component $!\n";
      }

      $linetitle = $s_string . $field_string;

    } elsif ($file =~ /s_([a-z]+)/) {

      if ($1 =~ /charge/) {
        $linetitle = q($g^{\pm}(k)$);
      } elsif ($1 =~ /direct/) {
        $linetitle = q($g^{\pm}(r)$);
      } else {
        die "File name doesn't match regex. $file $!\n";
      }

    } else {
      die "File name doesn't match anything. $file $!\n";
    }

    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
  }
} else {
  push @titles, '';
}

# works better for canonical
# my $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, charge value = $parameters{charge_value} * 2 {/Symbol p}.\n\n$meas measurements from $parameters{measurement_sweeps} MC steps.);


my $basedir = File::Spec->curdir();
my $plotpath = "$basedir/scripts/";
my $plotsuffix = ".tex";
my $gnuplotscript = "$plotpath" . "heatmap_latex.gp";

my $inpath = "$dir/";
my $insuffix = ".dat";
my $outpath = "$dir/plots";
my $tempsuffix = '.temp';

for my $i (0..$#filenames) {
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @tempfiles, $inpath . $filenames[$i] . $tempsuffix;
  push @outputfiles, $outpath . '/' . $filenames[$i] . $plotsuffix;

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

    # push @gnuplotargs, qq(FILE='$inputfiles[$i]'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$plottitle'; LINETITLE = '$titles[$i]'; PALETTE = '$palette';);
    push @gnuplotargs, qq(FILE='$inputfiles[$i]'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$titles[$i]'; LINETITLE = ''; PALETTE = '$palette';);
    push @latexargs, qq(latex -interaction=batchmode -output-directory=$outpath $outputfiles[$i]);
    push @dvipsargs, qq(dvips -q -D10000 -o $outpath/$filenames[$i].ps $outpath/$filenames[$i].dvi);
    push @ps2pdfargs, qq(ps2pdf -dPDFSETTINGS=/prepress -dColorImageResolution=600 $outpath/$filenames[$i].ps $outpath/$filenames[$i].pdf);

    #if (index($inputfiles[$i],'s_direct') != -1) {
    #  push @gnuplotargs, qq( PITICS = 'N');
    #} else {
    #  push @gnuplotargs, qq( PITICS = 'Y');
    #}

  }

  $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  system($syscall);
  system($latexargs[$i]);
  system($dvipsargs[$i]);
  system($ps2pdfargs[$i]);

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

sub commify {
  local $_  = shift;
  1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
  return $_;
}
