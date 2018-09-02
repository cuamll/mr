# set terminal cairolatex standalone size 15cm,15cm 
set terminal cairolatex standalone size 15cm,15cm header "\\usepackage[T1]{fontenc}\n\\renewcommand*\\familydefault{\\sfdefault}"
set output OUTPUT
set grid xtics ytics dt 2 lt 2 lc rgb '#666666'
set print "-"
set autoscale fix
load PALETTE

if (exists("PLOTTITLE")) {
  set title PLOTTITLE
}
set key top right spacing 3 width 2
set tmargin at screen 0.12
set bmargin at screen 0.88
set lmargin at screen 0.10
set rmargin at screen 0.86

if (exists("PITICS")) {
  set xtics pi;
  set ytics pi;
  set format x '%.0P';
  set format y '%.0P';
  set xlabel offset 0,-2 "$ G_x $";
  set ylabel offset -2,0 "$ G_y $";
  set format cb '$%g$'
}

# We have to add 0 to COLUMN, otherwise it's read as a string
COLUMN = COLUMN + 0

if (COLUMN==7) plot FILE u 1:2:($3 + $6) w image title LINETITLE; else \
  plot FILE u ($1/2):($2/2):COLUMN w image title LINETITLE;
