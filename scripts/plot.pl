#!/opt/local/bin/perl
# plot.pl: gnuplot plotting of output from maggs-rossetto CG code.
# called from analyse script with parameters included automatically
#
# Basically this script takes all the output files from the fortran code
# and generates a set of parameters which are then passed to gnuplot,
# and from there to latex -> dvips -> ps2pdf, in order to generate a set of
# labelled PDFs automatically. A lot of the script was taken up generating
# correct titles which render properly in latex, which makes it a lot less
# readable; I then turned off the titles to put the figures in my thesis,
# so this script could probably be tidied up a lot. But as it is I still
# have the option of turning the titles back on, so I've left it as is.

use strict;
use warnings;
use Env;
use Cwd;
use Getopt::Long;
use File::Path qw(make_path);
use File::Copy;
use File::Basename;

my $three_d = 0;
my $dir; my $kz; my @columns;
my $slots = 1;
my $palette = '~/.config/gnuplot/inferno.pal';
my $addtitles = 0; my $keep_aux = 0;
my @inputfiles; my @outputfiles; my @tempfiles;
my @gnuplotargs; my @latexargs; my @dvipsargs; my @ps2pdfargs;
my $input = GetOptions ("d=s"=> \$dir,
                     "t=i"=> \$addtitles,
                     "s=i"=> \$slots,
                     "k=i"=> \$keep_aux,
                     "p=s"=> \$palette);

if ($palette !~ /.*\.pal/) {
  $palette = "$palette.pal";
}

# we pass the absolute path to the timestamp directory. now to get parameters:
my %parameters = get_parameters("$dir/temp");

my $basedir = File::Spec->curdir();
my $plotpath = "$basedir/scripts/";
my $plotsuffix = ".tex";
my $gnuplotscript = "$plotpath" . "heatmap_latex.gp";

my $inpath = "$dir/";
my $insuffix = ".dat";
my $outpath = "$dir/plots";
print "Creating directory $outpath .\n";
make_path($outpath);
my $tempsuffix = '.temp';

my $steps = $slots * $parameters{no_samples} * $parameters{measurement_sweeps};
my $meas = $steps / $parameters{sample_interval};
my $steps_c = commify($steps);
my $meas_c = commify($meas);
my $chgen = lc $parameters{charge_generation};

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
my $plottitle_base;
if ($parameters{canon} =~ /T/ || $parameters{canon} =~ /Y/) {
  $plottitle_base = qq(Canonical: L = $parameters{L}, T = $parameters{temperature}, $parameters{charges} charges ($chgen), );
} else {
  $plottitle_base = qq(Grand canonical: L = $parameters{L}, T = $parameters{temperature}, \$ \\epsilon_c = $parameters{e_c} \$, );
}

# the s_ab_whatever files have four components; we want to plot each separately
for my $i (0..$#filenames) {
  my $file = $filenames[$i];
  if ($file =~ /s_ab_([a-z]+)/) {
    my $field_component = $1;

    # the multiplot is really awkward. ugly code. whatever
    my $mplotscript = "$plotpath" . "heatplot_multiplot.gnu";
    my $multiplottitle = qq(Grand canonical: L = $parameters{L}, T = $parameters{temperature}, \$ \\epsilon_c = $parameters{e_c} \$, );
    $linetitle = qq(\$ S^{\\alpha \\beta}_{$field_component} \$ );
    $multiplottitle = $multiplottitle . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    my $minfile = $inpath . $filenames[$i] . $insuffix;
    my $moutfile = $outpath . '/' . "s_ab_$field_component" . '_gnu';
    my $mplotarg = qq(FILE='$minfile'; OUTPUT='$moutfile$plotsuffix'; PALETTE = '$palette'; PITICS = 'Y';);
    if ($addtitles) {
      $mplotarg .= qq( PLOTTITLE='$multiplottitle';);
    }
    my $mlatexarg = qq(latex -interaction=batchmode -output-directory=$outpath $moutfile$plotsuffix > /dev/null);
    my $mdvipsarg = qq(dvips -q -D10000 -o $moutfile.ps $moutfile.dvi);
    my $mps2pdfarg =  qq(ps2pdf -dPDFSETTINGS=/prepress -dColorImageResolution=600 $moutfile.ps $moutfile.pdf 2> /dev/null);
    my $syscall = qq(gnuplot -e "$mplotarg" $mplotscript);
    print "Plotting $moutfile\n";
    system($syscall);
    system($mlatexarg);
    system($mdvipsarg);
    system($mps2pdfarg);
    if (!$keep_aux) {
      unlink("$moutfile.log");
      unlink("$moutfile.aux");
      unlink("$moutfile.dvi");
      unlink("$moutfile.ps");
      unlink("$moutfile-inc.eps");
    }

    # do each tensor component separately
    $linetitle = qq(\$ S^{xx}_{$field_component} \$ );
    $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_xx_$field_component";
    push @columns, 3;

    $linetitle = qq(\$ S^{xy}_{$field_component} \$ );
    $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_xy_$field_component";
    push @columns, 4;

    $linetitle = qq(\$ S^{yx}_{$field_component} \$ );
    $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_yx_$field_component";
    push @columns, 5;

    $linetitle = qq(\$ S^{yy}_{$field_component} \$ );
    $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
    push @titles, $plottitle;
    push @inputfiles, $inpath . $filenames[$i] . $insuffix;
    push @outputfiles, $outpath . '/' . "s_yy_$field_component";
    push @columns, 6;

    $linetitle = qq(\$ S^{xx}_{$field_component} + S^{yy}_{$field_component} \$ );
    $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
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
      print "File name doesn't match regex. $file $!\n";
      next;
    }

  } else {
    print "File name doesn't match anything. $file $!\n";
    next;
  }

  # if they're going in a document as a figure, might not want the titles
  $plottitle = $plottitle_base . qq($linetitle\n\n$meas_c measurements from $steps_c MC steps.);
  push @titles, $plottitle;

  # generate lists of input/output files and command-line arguments
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @outputfiles, $outpath . '/' . $filenames[$i];
  push @columns, 3;
}

warn "Different number of titles and files!\n" unless @inputfiles == @outputfiles;

for my $i (0..$#inputfiles) {
  # build up the command to run gnuplot with parameters
  # that match the data file we're plotting
  push @gnuplotargs, qq(FILE='$inputfiles[$i]'; OUTPUT='$outputfiles[$i]$plotsuffix'; COLUMN='$columns[$i]'; LINETITLE = ''; PALETTE = '$palette';);

  # now we need to run latex to get a DVI, dvips for a PS, ps2pdf for a PDF
  # so generate the correct set of commands to do this for a given output file
  push @latexargs, qq(latex -interaction=batchmode -output-directory=$outpath $outputfiles[$i]$plotsuffix > /dev/null);
  push @dvipsargs, qq(dvips -q -D10000 -o $outputfiles[$i].ps $outputfiles[$i].dvi);
  push @ps2pdfargs, qq(ps2pdf -dPDFSETTINGS=/prepress -dColorImageResolution=600 $outputfiles[$i].ps $outputfiles[$i].pdf 2> /dev/null);

  if (index($inputfiles[$i],'s_direct') != -1) {

  } else {
    $gnuplotargs[$i] .= qq( PITICS = 'Y';);
  }

  if ($addtitles) {
    $gnuplotargs[$i] .= qq( PLOTTITLE = '$titles[$i]';);
  }

  # because of the syntax of the gnuplot commands and the extra
  # switches just above, need to finally concatenate everything
  my $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  print "Plotting $outputfiles[$i] with column $columns[$i]\n";

  # now actually run gnuplot, latex, dvips, ps2pdf
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
