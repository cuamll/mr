set terminal pngcairo dashed enhanced font 'Helvetica, 14' size 800,800
set output OUTPUT
# set log y 2
# set grid ytics dt 2 lt 2 lc rgb '#666666'
# set border 3 back lc rgb '#666666' lt 2
set ytics ('0.125' 0.125, '0.25' 0.25, '0.5' 0.5)
set xtics (6, 12, 24, 36, 48, 60)
# set tics nomirror
set yrange [0.04:0.6]
set title PLOTTITLE
set key top right
set ylabel "avg. population"
set xlabel "number of charges"
set tmargin at screen 0.15
set bmargin at screen 0.85
set lmargin at screen 0.15
set rmargin at screen 0.85

set arrow from 6,0.125 to 60,0.125 nohead lc rgb '#C4145E' lw 2 dt 1 
set arrow from 6,0.125 to 60,0.125 nohead lc rgb '#A1D490' lw 2 dt 2 
set arrow from 6,0.250 to 60,0.250 nohead lc rgb '#44A4C2' lw 2 dt 2 
set arrow from 6,0.500 to 60,0.500 nohead lc rgb '#F2D33C' lw 2 dt 2 
set label "T_{I}, T_{IV}" at first 60, first 0.125
set label "T_{II}" at first 60, first 0.25
set label "T_{III}" at first 60, first 0.5

#/* the linestyles here are to make the colours match bjorgvin's paper */
plot FILE index 0 w lp ls 2 ps 2 title "T_{I}", \
       '' index 1 w lp ls 3 ps 2 title "T_{II}", \
       '' index 2 w lp ls 4 ps 2 title "T_{III}", \
       '' index 3 w lp ls 1 ps 2 title "T_{IV}"
