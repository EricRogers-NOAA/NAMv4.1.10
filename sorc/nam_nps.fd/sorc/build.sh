#!/bin/bash
cd NMMB_init/NPS
./clean

module use -a  ../../../modulefiles/wcossdell/
module load  v4.1.0_build

./conf dell 

module list

cp configure.nps_wcoss_dmpara configure.nps

OUTDIR_DELL=../../../../../exec/dell.exec/
OUTDIR=../../exec

mkdir -p $OUTDIR

./compile geogrid
cp -L geogrid.exe ${OUTDIR}/nam_geogrid

./compile metgrid
cp -L metgrid.exe ${OUTDIR}/nam_metgrid

./compile ungrib
cp -L ungrib.exe  ${OUTDIR}/nam_ungrib

./compile nemsinterp
cp -L nemsinterp.exe ${OUTDIR}/nam_nemsinterp

# this clean might not be necessary, but seems like a good idea
./clean

cp configure.nps_wcoss_serial configure.nps

./compile metgrid
cp -L metgrid.exe ${OUTDIR}/nam_metgrid_boco

./compile ungrib
cp -L ungrib/src/ungribs.exe  ${OUTDIR}/nam_ungrib_boco

# not needed
##./compile nemsinterp
##cp -L nemsinterp.exe ${OUTDIR}/nam_nemsinterp_boco

pwd
cp ${OUTDIR}/* ${OUTDIR_DELL}/.
cd ../..
