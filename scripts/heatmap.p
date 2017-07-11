set terminal pngcairo dashed enhanced font 'Helvetica, 14' size 800,800
set output OUTPUT
set view map
load PALETTE
set key center bottom
set grid xtics mxtics ytics dt 2 lt 2 lc rgb '#666666'
set xtics pi
set ytics pi
set format x '%.0Pπ'
set format y '%.0Pπ'
set tmargin at screen 0.10
set bmargin at screen 0.85
set lmargin at screen 0.15
set rmargin at screen 0.85
set title PLOTTITLE offset 0,2
splot FILE w image title LINETITLE
