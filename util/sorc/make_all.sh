#!/bin/bash
export MACHID=dell
echo 'Running on DELL'

if [ "${MACHID}" = dell ]; then
. /usrx/local/prod/lmod/lmod/init/bash
elif [ "${MACHID}" = theia ]; then
module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
fi
module purge
module use -a ../../modulefiles/${MACHID}
module load build/v4.0.0_build

export OUTmain=`dirname $(readlink -f ../exec/${MACHID}.exec/ )`
export OUTDIR=${OUTmain}/${MACHID}.exec

make -f ./Makefile

cd wgrib2_v2.0.7a/sorc
./install_all_grib_util_wcoss.sh_wgrib2only
cd ../exec
cp wgrib2 ../../../exec/dell.exec/.

set +x
module unload build/v4.0.0_build
set -x
exit 0
