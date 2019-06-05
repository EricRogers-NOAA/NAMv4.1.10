#!/bin/sh

########################################
# Script to clean up NAM model run directories 
########################################

set -xa

cd $DATA 

if [ ! -s ${GESDIR}/cleaned_${CDATE}_${RUNTYPE}.out ]; then

# Extract hourly precipitation from 00z NAM
# 12-36 h forecast for OCONUS soil moisture adjustment
# during NDAS
#
  #Workflow manager ought to ensure that the post/prdgen jobs are successful and so testing for
  # completion is not needed.
  # For NCO : these were moved to the exnam_prdgen script
# if [ $cyc -eq 00 ]; then
#   export CYCLE=${PDY}${cyc}
#   ${USHnam}/nam_getpcp_nam00z.sh
#   ${USHnam}/nam_getpcp_nam00z_nest.sh conus
# fi

  set -x

    if [ -d ${DATAROOT} ]; then
      cd ${DATAROOT}
      rm -rf ${DATAROOT}/${RUN}_${cyc}_pmgr_conus_${tmmark}_${envir}
      if [ $tmmark = tm00 ]; then
        #And clean up any remaining nam_${cyc}_tm*_${envir} directories that may have been missed previously
        rm -rf ${DATAROOT}/${RUN}_${cyc}_tm*_${envir}
        rm -rf ${DATAROOT}/${RUN}_${cyc}_pmgr_conus_tm*_${envir}
        rm -rf ${DATAROOT}/${RUN}_${cyc}_pmgr_alaska_tm*_${envir}
        rm -rf ${DATAROOT}/${RUN}_${cyc}_main_${envir}
      fi
      date
      cd $DATAROOT
    fi
    echo "ALL_CLEAN!" > ${GESDIR}/cleaned_${CDATE}_${RUNTYPE}.out
fi

exit 0 

