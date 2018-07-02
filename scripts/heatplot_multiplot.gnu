# #!/usr/bin/gnuplot
#
# set terminal cairolatex standalone size 15cm,15cm 
set terminal cairolatex standalone size 20cm,20cm header "\\usepackage[T1]{fontenc}\n\\renewcommand*\\familydefault{\\sfdefault}"
set output OUTPUT
set macros
set grid xtics ytics dt 2 lt 2 lc rgb '#666666'
set print "-"
set autoscale fix
load PALETTE

# set title PLOTTITLE
set key top right spacing 3 width 2
set tmargin at screen 0.12
set bmargin at screen 0.88
set lmargin at screen 0.12
set rmargin at screen 0.88


# x- and ytics for each row resp. column
NOXTICS = "unset xtics; \
           unset xlabel"

XTICS = "set xtics pi; \
         set format x '%.0P'; \
         set xlabel offset 0,-2 '$ G_x $'"

NOYTICS = "unset ytics; \
           unset ylabel"

YTICS = "set ytics pi; \
         set format y '%.0P'; \
         set ylabel offset -2,0 '$ G_y $'"

# Margins for each row resp. column
# TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.55"
# BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"
# LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55"
# RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95"
TMARGIN = "set tmargin at screen 0.88; set bmargin at screen 0.58"
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.25"
LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.40"
RMARGIN = "set lmargin at screen 0.48; set rmargin at screen 0.78"
CBTICS = "set cbtics scale 0.8; set format cb '%.2f'"

if (exists("PITICS")) {
  set xtics pi;
  set ytics pi;
  set format x '%.0P';
  set format y '%.0P';
  set xlabel offset 0,-2 "$ G_x $";
  set ylabel offset -2,0 "$ G_y $";
}

# We have to add 0 to COLUMN, otherwise it's read as a string
# COLUMN = COLUMN + 0

# if (COLUMN==7) plot FILE u 1:2:($3 + $6) w image title LINETITLE; else \
#   plot FILE u ($1/2):($2/2):COLUMN w image title LINETITLE;

set label PLOTTITLE at screen 0.5, 0.93 center front
@CBTICS

set multiplot layout 2,2 rowsfirst title ''
# --- GRAPH a
@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
plot FILE u ($1/2):($2/2):3 w image title '';
# --- GRAPH b
@TMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
plot FILE u ($1/2):($2/2):4 w image title '';
# --- GRAPH c
@BMARGIN; @LMARGIN
@XTICS; @YTICS
plot FILE u ($1/2):($2/2):5 w image title '';
# --- GRAPH d
@BMARGIN; @RMARGIN
@XTICS; @NOYTICS
plot FILE u ($1/2):($2/2):6 w image title '';
unset multiplot
### End multiplot
