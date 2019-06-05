#!/bin/sh

#------------------------------------------------------------
# Build all "emcsfc" programs.
#
# For more details, see the documentation in each
# program sub-directory.
#------------------------------------------------------------

#set -x

for directory in ./*.fd
do
  case $directory in
    *gridgen_sfc.fd)
      cd $directory/lib
      make clean
      make.sh
      cd ../driver
      make clean
      make.sh
      cd ../.. ;;
    *)
      cd $directory
      make clean
      make.sh
      cd .. ;;
  esac

done

cp ../exec/*  ../../../exec/dell.exec/.

echo; echo DONE BUILDING EMCSFC PROGRAMS
