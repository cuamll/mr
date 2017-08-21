# set term postscript eps size 18cm,18cm color colortext 12
set terminal cairolatex standalone size 20cm,20cm
# set terminal pngcairo dashed enhanced font 'Helvetica, 14' size 800,800
set output OUTPUT
# set log y 2
# set log x
# set grid ytics dt 2 lt 2 lc rgb '#666666'
# set border 3 back lc rgb '#666666' lt 2
# set xtics (0.25, 0.5, 1.0, 5.0)
# set tics nomirror
# set ytics ('0.125' 0.125, '0.25' 0.25, '0.5' 0.5)

XLOWER = 0.7
XUPPER = 2.0
YLOWER = 0
YUPPER = 2.0
POINTSIZE = 0.5

set yrange [YLOWER:YUPPER]
set xrange [XLOWER:XUPPER]
set title PLOTTITLE
set key top left spacing 3 width -2
set ylabel '\Large $ \frac{N}{k_B T} \left( \left< \bar{\mathbf{E}}^2 \right> - \left< \bar{\mathbf{E}} \right>^2 \right) $'
set xlabel '\Large Temperature'
set tmargin at screen 0.15
set bmargin at screen 0.85
set lmargin at screen 0.15
set rmargin at screen 0.85

#set arrow from XLOWER,0.125 to XUPPER,0.125 nohead lc rgb '#C4145E' lw 2 dt 1 
#set arrow from XLOWER,0.125 to XUPPER,0.125 nohead lc rgb '#A1D490' lw 2 dt 2 
#set arrow from XLOWER,0.250 to XUPPER,0.250 nohead lc rgb '#44A4C2' lw 2 dt 2 
#set arrow from XLOWER,0.500 to XUPPER,0.500 nohead lc rgb '#F2D33C' lw 2 dt 2 
#set label '$T_{I}$, $T_{IV}$' at first XUPPER + XUPPER / 50, first 0.125
#set label '$T_{II}$' at first XUPPER + XUPPER / 50, first 0.25
#set label '$T_{III}$' at first XUPPER + XUPPER / 25, first 0.5

#/* the linestyles here are to make the colours match bjorgvin's paper */
#plot FILE index 0 w lp ls 2 ps POINTSIZE title '\footnotesize L = 16', \
#       '' index 1 w lp ls 3 ps POINTSIZE title '\footnotesize L = 20', \
#       '' index 2 w lp ls 4 ps POINTSIZE title '\footnotesize 8 chg.', \
#       '' index 3 w lp ls 1 ps POINTSIZE title '\footnotesize 16 chg.', \
#       '' index 4 w lp ls 5 ps POINTSIZE title '\footnotesize 32 chg.'
plot FILE index 1 w lp ls 2 ps POINTSIZE title '\footnotesize L = 20'

#plot FILE using 1:4 index 0 w lp ls 2 ps POINTSIZE title '\footnotesize $T_{I}$', \
#       '' using 1:4 index 1 w lp ls 3 ps POINTSIZE title '\footnotesize $T_{II}$', \
#       '' using 1:4 index 2 w lp ls 4 ps POINTSIZE title '\footnotesize $T_{III}$', \
#       '' using 1:4 index 3 w lp ls 1 ps POINTSIZE title '\footnotesize $T_{IV}$', \
#       '' using 1:2 index 4 w lp ls 5 ps POINTSIZE title '\footnotesize $\left< \left| \bar{\mathbf{E}}^2 \right| \right> - \left< \left| \bar{\mathbf{E}} \right| \right>^2$', \
#       '' using 1:2 index 5 w lp ls 6 ps POINTSIZE title '\footnotesize Charge hop acceptance rate'
