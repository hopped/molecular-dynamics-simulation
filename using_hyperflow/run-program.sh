#!/bin/sh
# A script to compile and run the molecular docking program

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
        #cho "** found UUID = $UUID"
    else
        args="$args $var"
    fi
done  
#echo "UUID = $UUID -- args = $args"

DIR="/paasage"
cd $DIR/$UUID/Molecular_Docking

# NOTE: debugging
#echo
#echo $DIR/$UUID/Molecular_Docking
#ls -lah $DIR/$UUID/Molecular_Docking
#echo

sh compile.sh
sh run-cmd.sh $args  

