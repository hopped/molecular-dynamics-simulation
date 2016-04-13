#!/bin/sh
# Copyright 2013 University of Stuttgart, Germany
# Author: Anthony Sulistio (HLRS)
#
# A script to compile the source code
# usage: ./compile.sh [debug]

VAL="RELEASE"
if [ -n "$1" ]; then
    VAL="DEBUG"
fi

cd src
make distclean 
make ARCH=CPU TARGET=$VAL MPI=1 COMPILER=mpi -j4
cp CMD_CPU main

