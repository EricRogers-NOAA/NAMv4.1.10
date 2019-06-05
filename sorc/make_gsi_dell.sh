#!/bin/bash

set -ex

cd ..
pwd=$(pwd)
dir_root=$pwd

. $MODULESHOME/init/sh
conf_target=nco

[ -d $dir_root/exec ] || mkdir -p $dir_root/exec

dir_modules=$dir_root/modulefiles
module purge
module use $dir_modules/
module load modulefile.nam_gsi.wcoss_d
module list

dir_refsorc=$dir_root/sorc/nam_ref2nemsio.fd
cd $dir_refsorc
export OUTDIR=$dir_root/exec/dell.exec
make -f makefile

dir_libsorc=$dir_root/sorc/nam_gsi.fd/lib/GSD/gsdcloud4nmmb
cd $dir_libsorc
make -f makefile clean
make -f makefile

dir_sorc=$dir_root/sorc/nam_gsi.fd/src
cd $dir_sorc
make -f Makefile clean
make -f Makefile -j 4
cp -p nam_gsi $dir_root/exec/dell.exec

