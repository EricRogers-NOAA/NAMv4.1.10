#!/bin/bash

SORCnam=$(pwd)
echo 'Running on DELL'
MACHID=dell

if [ "${MACHID}" = dell ]; then
. /usrx/local/prod/lmod/lmod/init/bash
fi
module purge
module use -a ../modulefiles/${MACHID}
module load build/v4.0.0_build

set -x
export OUTmain=`dirname $(readlink -f ../exec/${MACHID}.exec/ )`
export OUTDIR=${OUTmain}/${MACHID}.exec

make -f ./Makefile

# Now build the GSI 
./make_gsi_dell.sh > build_gsi.log 2>> build_gsi.log

# Now build ncep post
cd ${SORCnam}/nam_nceppost.fd/sorc
./build_ncep_post.sh > build_post.log 2>> build_post.log

# Now build NMMB forecast model code
cd ${SORCnam}/nam_nems_nmmb.fd/src
./build.sh > build_fcst.log 2>> build_fcst.log

# Now build land-surface utility codes
cd ${SORCnam}/nam_land_utilities.fd/sorc
./build_emcsfc.sh > build_emcsfc.log 2>> build_emcsfc.log

cd $SORCnam

echo -e " I do not build the following external source code: \n" \
     "   nam_nps.fd "
echo " Please compile this using the ./nam.v4.1.10/sorc/nam_nps.fd/sorc/build.sh script."

exit 0
