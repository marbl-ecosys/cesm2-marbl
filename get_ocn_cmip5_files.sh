#!/bin/bash

# Scipt by klindsay, this is adaptation of
# https://github.com/klindsay28/CESM2_coup_carb_cycle_JAMES/blob/master/get_atm_cmip5_files.sh

if [ $# -ne 1 ]; then
  echo "ERROR: this script requires a single argument"
  echo "Usage: $0 VAR_NAME"
  exit 1
fi

VAR_NAME=$1

#for case in b40.prescribed_carb.001 \
#            b40.1850_ramp.1deg.ncbdrd.001 \
#            b40.1850_ramp.1deg.ncbdrc.001 \
#            b40.1850_ramp.1deg.ncbcrd.002 \
#            b40.20th.1deg.bdrd.001 \
#            b40.coup_carb.004 \
#            b40.20th.1deg.coup.001 ; do
for case in b40.20th.1deg.bdrd.001 ; do

   echo $case
   cd /glade/p/cgd/oce/projects/cesm2-marbl/intake-esm-data/$case/ocn/proc/tseries/monthly
   pwd
   hsi "cd /CCSM/csm/$case/ocn/proc/tseries/monthly ; ls *${VAR_NAME}*"
   hsi "cd /CCSM/csm/$case/ocn/proc/tseries/monthly ; cget *${VAR_NAME}*"

done

