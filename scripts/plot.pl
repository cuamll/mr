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
my $dir; my $kz; my @columns;
my $palette = '~/.config/gnuplot/inferno.pal';
my $addtitles = 1; my $keep_aux = 0;
my @inputfiles; my @outputfiles; my @tempfiles;
my @gnuplotargs; my @latexargs; my @dvipsargs; my @ps2pdfargs;
my $input = GetOptions ("d=s"=> \$dir,
                     "t=i"=> \$addtitles,
                     "k=i"=> \$keep_aux,
                     "p=s"=> \$palette);

if ($palette !~ /.*\.pal/) {
  $palette = "$palette.pal";
}

# we pass the absolute path to the timestamp directory. now to get parameters:
my %parameters = get_parameters("$dir/input.in");

my $basedir = File::Spec->curdir();
my $plotpath = "$basedir/scripts/";
my $plotsuffix = ".tex";
my $gnuplotscript = "$plotpath" . "heatmap_latex.gp";

my $inpath = "$dir/";
my $insuffix = ".dat";
my $outpath = "$dir/plots";
make_path($outpath);
my $tempsuffix = '.temp';

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
my $s_string; my $field_string; my $linetitle; my $plottitle;

# the s_ab_whatever files have four components; we want to plot each separately
for my $i (0..$#filenames) {
  my $file = $filenames[$i];
  if ($file =~ /s_ab_([a-z]+)/) {
    my $field_component = $1;

    # do each tensor component separately
    $linetitle = qq(\$ S^{xx}_{$field_component} \$ );
    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_xx_$field_component";
    push @columns, 3;

    $linetitle = qq(\$ S^{xy}_{$field_component} \$ );
    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_xy_$field_component";
    push @columns, 4;

    $linetitle = qq(\$ S^{yx}_{$field_component} \$ );
    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_yx_$field_component";
    push @columns, 5;

    $linetitle = qq(\$ S^{yy}_{$field_component} \$ );
    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_yy_$field_component";
    push @columns, 6;

    $linetitle = qq(\$ S^{xx}_{$field_component} + S^{yy}_{$field_component} \$ );
    $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_trace_$field_component";
    push @columns, 7;
  }
}

for my $i (0..$#filenames) {
  my $file = $filenames[$i];

  # already done these
  if ($file =~ /s_ab_([a-z]+)/) {
    next;
  }

  if ($file =~ /s_([a-z]+)_([a-z]+)/) {
    my $s_component = $1;
    my $field_component = $2;

    if ($s_component =~ /xx/) {
      $s_string = qq(\$ S^{xx}_{$field_component} \$);
    } elsif ($s_component =~ /perp/) {
      $s_string = qq(\$ S^{\\perp}_{$field_component} \$);
    } elsif ($s_component =~ /par/) {
      $s_string = qq(\$ S^{\\parallel}_{$field_component} \$);
    } else {
      die "s_component is wrong: $s_component $!\n";
    }

    $linetitle = $s_string;

  } elsif ($file =~ /s_([a-z]+)/) {

    if ($1 =~ /charge/) {
      $linetitle = q($ g^{\pm}(k) $);
    } elsif ($1 =~ /direct/) {
      $linetitle = q($ g^{\pm}(r) $);
    } else {
      die "File name doesn't match regex. $file $!\n";
    }

  } else {
    die "File name doesn't match anything. $file $!\n";
  }

  # if they're going in a document as a figure, might not want the titles
  $plottitle = qq(L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges, $linetitle\n\n$meas_c measurements from $steps_c MC steps.);
  push @titles, $plottitle;

  # generate lists of input/output files and command-line arguments
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @outputfiles, $outpath . '/' . $filenames[$i];
  push @columns, 3;
}

warn "Different number of titles and files!\n" unless @inputfiles == @outputfiles;

for my $i (0..$#inputfiles) {
  push @gnuplotargs, qq(FILE='$inputfiles[$i]'; OUTPUT='$outputfiles[$i]$plotsuffix'; COLUMN='$columns[$i]'; LINETITLE = ''; PALETTE = '$palette';);
  push @latexargs, qq(latex -interaction=batchmode -output-directory=$outpath $outputfiles[$i]$plotsuffix > /dev/null);
  push @dvipsargs, qq(dvips -q -D10000 -o $outputfiles[$i].ps $outputfiles[$i].dvi);
  push @ps2pdfargs, qq(ps2pdf -dPDFSETTINGS=/prepress -dColorImageResolution=600 $outputfiles[$i].ps $outputfiles[$i].pdf);

  if (index($inputfiles[$i],'s_direct') != -1) {
    # $gnuplotargs[$i] .= qq( PITICS = 'N';);
  } else {
    $gnuplotargs[$i] .= qq( PITICS = 'Y';);
  }

  if ($addtitles) {
    $gnuplotargs[$i] .= qq( PLOTTITLE = '$titles[$i]';);
  }

  my $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  print "Plotting $outputfiles[$i] with column $columns[$i]\n";
  system($syscall);
  system($latexargs[$i]);
  system($dvipsargs[$i]);
  system($ps2pdfargs[$i]);

  if (!$keep_aux) {
    unlink("$outputfiles[$i].log");
    unlink("$outputfiles[$i].aux");
    unlink("$outputfiles[$i].dvi");
    unlink("$outputfiles[$i].ps");
    unlink("$outputfiles[$i]-inc.eps");
  }

}

# warn "Different number of titles and files!\n" unless @titles == @filenames;

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
