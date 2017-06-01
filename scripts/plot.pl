#!/opt/local/bin/perl
# PRELIMINARY!
# i'm probably gonna expand the scripting to automate more stuff.
# currently this just plots heatmaps of correlation functions
# call with "perl plot.pl -l=L ..." etc. from scripts dir
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $length;
my $meas;
my $steps;
my $chg;
my $kz;
my $fileprefix;
my $outputpath;
my @inputfiles;
my @outputfiles;
my @tempfiles;
my @gnuplotargs;
my $input;
$input = GetOptions ("l=i"=> \$length,
                     "m=i"=> \$meas,
                     "s=s"=> \$steps,
                     "c=i"=> \$chg,
                     "k=i"=> \$kz,
                     "fp=s"=> \$fileprefix,
                     "o=s"=> \$outputpath);
my $plottitle = "L = $length, $meas measurements from " .
                "$steps MC steps, $chg charges";

# construct the input and output files to pass to gnuplot
# call this script from $MR_DIR! so curdir() gives the base directory
my @filenames = ('charge_struc',
                 'alt_charge_struc',
                 'field_struc',
                 'alt_field_struc',
                 'irrot_field_struc',
                 's_perp',
                 'alt_s_perp_total');

my $basedir = File::Spec->curdir();
my $inpath = "$basedir/out/";
my $insuffix = '.dat';
my $tempsuffix = '.temp';
my $plotpath = "$basedir/plots/";
my $plotsuffix = ".png";
my $gnuplotscript = "$plotpath" . "heatmap.p";

# get the relevant lines of the plot files based on kz
my $lowerbound = (((2 * $length) + 1)**2) * ($kz + $length);
my $upperbound = $lowerbound + 1 + ((2 * $length) + 1)**2;

for my $i (0..$#filenames) {
  # run awk to create temp files
  my $awkcall;
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @tempfiles, $inpath . $filenames[$i] . $tempsuffix;
  push @outputfiles, $outputpath . $filenames[$i] . $plotsuffix;

  $awkcall = qq[awk '{ if (NR>$lowerbound && NR<$upperbound) print \$1,\$2,\$4 }' $inputfiles[$i] > $tempfiles[$i]];
  system($awkcall);

}

my @linetitles = ("Charge-charge structure factor at k_z = $kz",
                  "Charge-charge structure factor (read from unformatted file) at k_z = $kz",
                  "Field-field structure factor at k_z = $kz",
                  "Field-field structure factor (read from unformatted file) at k_z = $kz",
                  "Field-field structure factor - irrotational - at k_z = $kz",
                  "S_{⟂} at k_z = $kz",
                  "S_{⟂} (read from unformatted file) at k_z = $kz");

warn "Different number of titles and files!\n" unless @linetitles == @filenames;

for my $i (0..$#filenames) {
  my $syscall;
  push @gnuplotargs, qq(FILE='$tempfiles[$i]'; KZ = '$kz'; LSIZE = '$length'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$plottitle'; LINETITLE = '$linetitles[$i]');

  $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  system($syscall);

  # delete the temp files
  unlink($tempfiles[$i]);
}
