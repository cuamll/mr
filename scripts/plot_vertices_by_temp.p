set terminal pngcairo dashed enhanced font 'Helvetica, 14' size 800,800
set output OUTPUT
# set log y 2
# set log x
# set grid ytics dt 2 lt 2 lc rgb '#666666'
# set border 3 back lc rgb '#666666' lt 2
# set xtics (0.25, 0.5, 1.0, 5.0)
# set tics nomirror
set ytics ('0.125' 0.125, '0.25' 0.25, '0.5' 0.5)

XLOWER = 0.05
XUPPER = 0.5
YLOWER = 0.04
YUPPER = 0.6

set yrange [YLOWER:YUPPER]
set xrange [XLOWER:XUPPER]
set title PLOTTITLE
set key top right
set ylabel "Avg. population"
set xlabel "Temperature"
set tmargin at screen 0.15
set bmargin at screen 0.85
set lmargin at screen 0.15
set rmargin at screen 0.85

set arrow from XLOWER,0.125 to XUPPER,0.125 nohead lc rgb '#C4145E' lw 2 dt 1 
set arrow from XLOWER,0.125 to XUPPER,0.125 nohead lc rgb '#A1D490' lw 2 dt 2 
set arrow from XLOWER,0.250 to XUPPER,0.250 nohead lc rgb '#44A4C2' lw 2 dt 2 
set arrow from XLOWER,0.500 to XUPPER,0.500 nohead lc rgb '#F2D33C' lw 2 dt 2 
set label "T_{I}, T_{IV}" at first XUPPER + XUPPER / 50, first 0.125
set label "T_{II}" at first XUPPER + XUPPER / 50, first 0.25
set label "T_{III}" at first XUPPER + XUPPER / 25, first 0.5

#/* the linestyles here are to make the colours match bjorgvin's paper */
plot FILE using 1:4 index 0 w lp ls 2 ps 2 title "T_{I}", \
       '' using 1:4 index 1 w lp ls 3 ps 2 title "T_{II}", \
       '' using 1:4 index 2 w lp ls 4 ps 2 title "T_{III}", \
       '' using 1:4 index 3 w lp ls 1 ps 2 title "T_{IV}"
