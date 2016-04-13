#!/bin/sh
# A script to set up webserver and NFS storage
# For setting up the webserver, the script needs to be run on the webserver itself
#
# For setting up NFS: setup_NFS.sh NFS_IP [start_IP] [end_IP] 
# Example: ./setup_NFS.sh 109.231.122.193  170  220 -> will create an NFS server 
# on 109.231.122.193 and scan IP from 109.231.122.170 to 109.231.122.220 for
# list of MPI worker nodes.

echo "Running $0" "$@"

# create UUID. Need uuid-runtime package in debian
UUID=`uuidgen`
args=""
i=0
foundID=false

# loop through each input parameters
array="$@"
for var in $array
do
    i=`expr $i + 1`
    #echo "i = $i -- $var"

    # get the directory name, which is based on the UUID
    if [ "$var" = "--dir" ]
    then
        foundID=true
        continue
    fi

    if [ $foundID = true ]; then
        UUID=$var
        foundID=false
        #echo "** found UUID = $UUID"
    else
        args="$args $var"
    fi
done  

## Setting up Webserver and NFS
#echo "UUID = $UUID -- args = $args"
#sh /home/paasage/script/setup_Webserver.sh
#sh /home/paasage/script/setup_NFS.sh $args

DIR="/paasage"
mkdir -vp $DIR/$UUID

cd $DIR/$UUID
cp -v /home/paasage/Molecular_Docking.tar.gz .
tar zxf Molecular_Docking.tar.gz
rm -f Molecular_Docking.tar.gz

