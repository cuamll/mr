# set terminal cairolatex standalone size 15cm,15cm 
set terminal cairolatex standalone size 14cm,12cm header "\\usepackage[T1]{fontenc}\n\\usepackage[fleqn]{amsmath}\n\\renewcommand*\\familydefault{\\sfdefault}"
set output OUTPUT
set grid xtics ytics dt 2 lt 2 lc rgb '#666666'
set print "-"
set autoscale fix
load PALETTE
# set cbrange [0:2.0]
# set logscale cb

if (exists("PLOTTITLE")) {
  set title PLOTTITLE
}
set key top right spacing 3 width 2
# plot is 8/7 times as wide as it is tall, to allow for colour bar
# so we need to make sure the screen ratio is 7/8 to make it square
set tmargin at screen 0.13
set bmargin at screen 0.93
set lmargin at screen 0.12
set rmargin at screen 0.82

if (exists("PITICS")) {
  set xtics pi offset 0,-0.5;
  set ytics pi offset -1.75,0;
  set format x '{\Large $ \mathsf{ %- .0P } $}';
  set format y '{\Large $ \mathsf{ %- .0P } $}';
  set xlabel offset 0,-1.5 "{\\Large $ \\mathsf{G_x} $}";
  set ylabel offset 3.5,0 "{\\Large $ \\mathsf{G_y} $}";
  set format cb '{\Large \raggedright $ \mathsf{ %g } $}'
  # set format cb '{\Large $ \mathsf{ % 5.3f } $}'
  set cbtics offset 3.5,0
}

# We have to add 0 to COLUMN, otherwise it's read as a string
COLUMN = COLUMN + 0

if (COLUMN==7) plot FILE u ($1/2):($2/2):($3 + $6) w image title LINETITLE; else \
  plot FILE u ($1/2):($2/2):COLUMN w image title LINETITLE;
