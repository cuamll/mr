#!/bin/bash

# there are definitely neater ways to do this
# e.g. change underscore to space in filename, use that

gnuplot -e "FILE='out/charge_struc.dat'; OUTPUT='plots/charge_struc.png'; PLOTTITLE='L = 6, 600 measurements from 15,000 MC steps, 54 charges'; LINETITLE='Charge-charge structure factor at k_z = 0" plots/heatmap.p

gnuplot -e "FILE='out/field_struc.dat'; OUTPUT='plots/field_struc.png'; PLOTTITLE='L = 6, 600 measurements from 15,000 MC steps, 54 charges'; LINETITLE='Field-field structure factor at k_z = 0" plots/heatmap.p

gnuplot -e "FILE='out/s_perp.dat'; OUTPUT='plots/s_perp.png'; PLOTTITLE='L = 6, 600 measurements from 15,000 MC steps, 54 charges'; LINETITLE='S_{⟂} at k_z = 0" plots/heatmap.p
