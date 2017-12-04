set terminal cairolatex standalone size 15cm,15cm
set output OUTPUT
set grid xtics ytics dt 2 lt 2 lc rgb '#666666'
set print "-"
load PALETTE

set title PLOTTITLE
set key top right spacing 3 width 2
set tmargin at screen 0.15
set bmargin at screen 0.85
set lmargin at screen 0.15
set rmargin at screen 0.85

if (exists("PITICS")) {
  set xtics pi;
  set ytics pi;
  set format x '%.0P';
  set format y '%.0P'
}

# We have to add 0 to COLUMN, otherwise it's read as a string
COLUMN = COLUMN + 0

if (COLUMN==7) plot FILE u 1:2:($3 + $6) w image title LINETITLE; else \
  plot FILE u 1:2:COLUMN w image title LINETITLE;
