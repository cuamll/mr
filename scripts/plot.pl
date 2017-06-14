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
                 'field_struc',
                 'irrot_field_struc',
                 'rot_field_struc',
                 's_perp',
                 'rot_s_perp',
                 'irrot_s_perp');

my @linetitles = ("Charge-charge structure factor at k_z = $kz",
                  "Field-field structure factor at k_z = $kz",
                  "Field-field structure factor - irrotational - at k_z = $kz",
                  "Field-field structure factor - rotational - at k_z = $kz",
                  "S_{⟂} at k_z = $kz",
                  "S_{⟂} - rotational - at k_z = $kz",
                  "S_{⟂} - irrotational - at k_z = $kz");

my $basedir = File::Spec->curdir();
my $inpath = "$basedir/out/";
my $insuffix = '.dat';
my $tempsuffix = '.temp';
my $plotpath = "$basedir/plots/";
my $plotsuffix = ".png";
my $gnuplotscript = "$plotpath" . "heatmap.p";

# get the relevant lines of the plot files based on kz
my $bz = 2;
my $sperp_size = 5;
my $lowerbound = ((($bz * $length) + 1)**2) * ($kz + ($bz * ($length / 2)));
my $upperbound = $lowerbound + 1 + (($bz * $length) + 1)**2;
my $splb = ((($sperp_size * $length) + 1)**2) * ($kz + ($sperp_size * ($length / 2)));
my $spub = $splb + 1 + (($sperp_size * $length) + 1)**2;

for my $i (0..$#filenames) {
  # run awk to create temp files
  my $awkcall;
  push @inputfiles, $inpath . $filenames[$i] . $insuffix;
  push @tempfiles, $inpath . $filenames[$i] . $tempsuffix;
  push @outputfiles, $outputpath . $filenames[$i] . $plotsuffix;

  if (index($inputfiles[$i],'s_perp') != -1) {
    $awkcall = qq[awk '{ if (NR>$splb && NR<$spub) print \$1,\$2,\$4 }' $inputfiles[$i] > $tempfiles[$i]];
  } else {
    $awkcall = qq[awk '{ if (NR>$lowerbound && NR<$upperbound) print \$1,\$2,\$4 }' $inputfiles[$i] > $tempfiles[$i]];
  }
  system($awkcall);

}

warn "Different number of titles and files!\n" unless @linetitles == @filenames;

for my $i (0..$#filenames) {
  my $syscall;
  push @gnuplotargs, qq(FILE='$tempfiles[$i]'; KZ = '$kz'; LSIZE = '$length'; OUTPUT='$outputfiles[$i]'; PLOTTITLE = '$plottitle'; LINETITLE = '$linetitles[$i]');

  $syscall = qq(gnuplot -e "$gnuplotargs[$i]" $gnuplotscript);
  system($syscall);

  # delete the temp files
  unlink($tempfiles[$i]);
}
