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
rm -f *.pov *.dat *.xyz

## Number of CPUs is automatically set via 'nproc'
CPU=$(nproc)

## number of molecules
NUM=1000
if [ -n "$1" ]; then
    NUM=$1
fi

## how long it will run, i.e. the time step is 0.000 0.005 0.010 0.015 ... end time
END="0.05"
if [ -n "$2" ]; then
    END=$2
fi

TEMPERATURE=0.85
if [ -n "$3" ]; then
    TEMPERATURE=$3
fi

OUTPUT_FILE=md-simulation.tgz
if [ -n "$4" ]; then
    OUTPUT_FILE=$4
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
./main -N $NUM -n 1 -T $TEMPERATURE --domain-type=cuboid --cutoff-radius=1 -m 1 \
--simulation-end-time=$END --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output

mpirun -np $CPU $INCL_HOST \
./main -N $NUM -n 1 -T $TEMPERATURE --domain-type=cuboid --cutoff-radius=1 -m 1 \
--simulation-end-time=$END --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output

tar zcvf $OUTPUT_FILE psp-*

mv $OUTPUT_FILE ..
