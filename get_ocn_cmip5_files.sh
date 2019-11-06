#!/bin/bash

for case in b40.prescribed_carb.001 \
            b40.1850_ramp.1deg.ncbdrd.001 \
            b40.1850_ramp.1deg.ncbdrc.001 \
            b40.1850_ramp.1deg.ncbcrd.002 \
            b40.20th.1deg.bdrd.001 \
            b40.coup_carb.004 \
            b40.20th.1deg.coup.001 ; do

   echo $case
   cd /glade/p/cgd/oce/projects/cesm2-marbl/intake-esm-data/$case/ocn/proc/tseries/monthly
   pwd
   hsi "cd /CCSM/csm/$case/ocn/proc/tseries/monthly ; ls *IRON_FLUX*"
   hsi "cd /CCSM/csm/$case/ocn/proc/tseries/monthly ; cget *IRON_FLUX*"
   # hsi "cd /CCSM/csm/$case/atm/proc/tseries/monthly ; cget *CO2* *SFCO2* *TMCO2* *TS.*"

done

