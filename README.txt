CMD

Copyright (y) 2011-2014 Christoph Niethammer <christoph.niethammer@gmail.com>

CMD is a Molecular Dynamics simulation program for HPC systems.
It is designed for high efficiency on CPUs as well as scalability on HPC systems.

License: the "New BSD License"

Requirements:

core:
gcc >= 4.6.1
glibc >= 2.14

optional:
* MPI library >= MPI-2: MPI parallelization
* libxml2: Suport for xml input files
* perl: make testrun


Directories:

doc:
examples:
src: source code
utils:


---------------------------------------------------------------------------
Package Dependencies

on Debian-based distro like Ubuntu:
apt-get install gcc make gfortran build-essential g++ openmpi-doc openmpi-bin libopenmpi-dev libc6-dev screen

Note: to install newer version of glibc on Debian 7
vim /etc/apt/sources.list
deb http://ftp.de.debian.org/debian testing main contrib non-free
apt-get -t testing install libc6-dev

---------------------------------------------------------------------------
For running on the 32-bit machine 

modify src/typedefs.h line 16 (uncomment the below line):
/* #define USE_FLOAT */     /** Note: running on 32-bit OS **/  

---------------------------------------------------------------------------
When running this MD application using Hyperflow (https://github.com/dice-cyfronet/hyperflow), 
here are the steps:
sudo mkdir /paasage/
sudo chmod 777 /paasage/

cd using_hyperflow/
vim Wf_Molecular_docking.json
-- change the location of the *.sh files in 
        "executable": "/home/paasage/script/...."

-- change the simulation parameters in
        "args": "2 1000 0.1"
   for  num_cores  num_molecules  end_simulation_time, respectively

vim pre-processing.sh
-- change the location of this MD application:
        cp -v /home/paasage/Molecular_Docking.tar.gz .

---------------------------------------------------------------------------
Compiling the CMD program

./compile.sh
    -- OR -- using MPI    
cd src
make distclean; make TARGET=RELEASE MPI=1 ARCH=CPU COMPILER=mpi -j4
cp CMD_CPU main 
    -- OR -- using OpenMP
cd src
export OMP_NUM_THREADS=4
make distclean; make TARGET=RELEASE OMP=1 ARCH=CPU COMPILER=gcc
cp CMD_CPU main 

---------------------------------------------------------------------------
Running the CMD program

To run using 4 cores and 1000 molecules using MPI:
./run-cmd.sh 4 1000 
    -- OR --
cd src
mpirun -np 4 ./main -N 1000 -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 --simulation-end-time=0.05 --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output


To run using 4 cores and 1000 molecules until 2 simulation time: 
./run-cmd.sh 4 1000 2
    -- OR --
cd src
mpirun -np 4 ./main -N 1000 -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 --simulation-end-time=2 --molecule-container=BASICN2 --thermostat=velocity-scaling --ascii-output --povray-output


========== For OpenMP only

./main -v -N 1000 -n 0.9 -T 0.85 --domain-type=cube --timestep-length=0.005 --cutoff-radius=2 -m 1 --simulation-end-time=1 --molecule-container=BASICN2 --thermostat=velocity-scaling --simulation-equilibration-time=0.5 --gridgenerator-lattice-centering=primitive --n-per-subdomain=20 --generator=dropOverBasin --generator-drop-radius=2 --ascii-output --povray-output


---------------------------------------------------------------------------

-n --rho : density
--verbose=5  : will generate the main.log file (value from 3 onwards)
--simulation-end-time=0.05 --> how long until it ends. The time step is 0.000 0.005 0.010 0.015 ... 0.05


For the list of parameters: ./main -h
main [[OPTIONS]...] [[FILE]...] 
   -N, --num-molecules=ULONG
       number of molecules
   -V, --volume=DOUBLE
   --domain-type=TYPE
         allowed types: cube, rectangle
   -T, --temperature=DOUBLE
   -n, --rho=DOUBLE
   --simulation-time=DOUBLE
   -v, --verbose=INT
   --cutoff-radius=DOUBLE
   --simulation-end-time=DOUBLE
   --simulation-start-time=DOUBLE
   --thermostat=TYPE
          allowed types: velocity-scaling
   --povray-output
   --ascii-output
   --molecule-container=TYPE
       allowed types: MOLECULEBLOCKS, BASICN2
   -h, --help

NOTE:
* If mpirun suddenly got killed using --ascii-output --povray-output, then
reduce the num of molecules or N to 600,000 or less
* To disable papi, comment all code related to PAPI in src/main.c file.
* To enable papi, uncomment src/main.c line 35: 
//#include <papi.h>

and in src/makefiles/mpi.mk link the papi lib manually
LDFLAGS_RELEASE = -O3 -L $(shell echo $(MPI_DIR))/lib -lmpi -L /usr/local/papi-3.7.2/lib -lpapi

------------------------------------------------------------------
Using POV-RAY for converting the results into PNG files

Installing povray:
wget http://www.povray.org/redirect/www.povray.org/ftp/pub/povray/Old-Versions/Official-3.62/Unix/povray-3.6.tar.gz
tar zxvf povray-3.6.tar.gz
apt-get install gawk build-essential autoconf automake
cd povray-3.6.1
./configure COMPILED_BY="your name <email@address>"
make
make install

To convert the POV file into PNG, type:  
./make-image.sh
    -- OR --
cd src    
cp -vf ../pov-template.inc psp-header.inc
WIDTH=1024  
HEIGHT=768
povray -w$WIDTH -h$HEIGHT +A -D +WL0 -GA file.pov

------------------------------------------------------------------
To create a movie from PNG files

apt-get install mencoder
cd src    
mencoder mf://*.png -mf w=$WIDTH:h=$HEIGHT:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.avi

------------------------------------------------------------------
To view the xyz file, on small molecules where N = 100 use: 
jmol avogadro gdis pymol xmakemol

http://jmol.sourceforge.net/
Note: gdis has a better GUI compared to others

