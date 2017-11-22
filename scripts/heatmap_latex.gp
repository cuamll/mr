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

#if (PITICS eq "Y") {
  set xtics pi
  set ytics pi
  set format x '%.0Pπ'
  set format y '%.0Pπ'
#} else {
#  set xtics
#  set ytics
#}

plot FILE w image title LINETITLE
