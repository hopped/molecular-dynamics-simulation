#!/bin/sh
# Copyright 2013 University of Stuttgart, Germany
# Author: Anthony Sulistio (HLRS)
#
# usage: ./run.sh [num_cpu] [num_molecule] [end_time]
#
# Examples:
# To run using 4 cores and 1000 molecules: ./run-cmd.sh 4 1000
# To run using 4 cores and 1000 molecules until 2 simulation time: 
# ./run-cmd.sh 4 1000 2

cp -vf pov-template.inc src/psp-header.inc
cd src/
rm -f *.pov *.dat *.xyz
#cp -vf ../pov-template.inc psp-header.inc

CPU=4
if [ -n "$1" ]; then
    CPU=$1
fi

## number of molecules
NUM=10000
if [ -n "$2" ]; then
    NUM=$2
fi

## how long it will run, i.e. the time step is 0.000 0.005 0.010 0.015 ... end time
END="0.05"
if [ -n "$3" ]; then
    END=$3
fi


## check if hostfile exists
INCL_HOST=""
HOSTFILE="hostfile.txt"
if [ -e $HOSTFILE ]; then
    INCL_HOST="-hostfile $HOSTFILE"
fi
echo "host file =" $INCL_HOST

# to compile
#make distclean; make ARCH=CPU TARGET=RELEASE MPI=1 COMPILER=mpi -j8
echo "Running the molecular docking experiment"

echo mpirun -np $CPU $INCL_HOST \
./main -N $NUM -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 \
--simulation-end-time=$END --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output

mpirun -np $CPU $INCL_HOST \
./main -N $NUM -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 \
--simulation-end-time=$END --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output

