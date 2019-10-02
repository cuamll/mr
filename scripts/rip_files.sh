#!/bin/bash

mrdir=/Users/cgray/code/mr/
for dir in ${mrdir}/out/salviati/*
do
  # obviously change this to do canonical results
  if [[ $dir =~ "gce" ]]
  then
    # ensures the cut data is actually there - 
    # can turn off if analysis is already done
    ./scripts/analyse.pl -spr=0 -quadrics=0 -helmholtz=1 -plot=0 -multiplot=0 -lorentz=0 -quiver=0 -d=$dir
    rho_t=$(grep -A 1 "Avg. charge" $dir/sphe_sus.dat)
    temp=$(printf "$rho_t" | sed -e '/#.*/d' | awk '{print $1}')
    rho=$(printf "$rho_t" | sed -e '/#.*/d' | awk '{print $2}')
    qdata=$(grep -A 129 "Simulated" $dir/helmholtz_twolor/fit_params.dat | sed -e 's/\(Sim.*\)/# \1/' -e 's/.*\[//' -e 's/\].*//')
    echo "T = $temp; printing to file T_${temp}_trace.dat"
    printf "$rho_t\n\n$qdata" > "T_${temp}_trace.dat"
  fi
done
