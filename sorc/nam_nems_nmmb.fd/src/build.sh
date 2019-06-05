#!/bin/bash

make clean
./configure dell
. conf/modules.nems.dell
./compile.sh_pll
