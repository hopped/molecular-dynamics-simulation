#!/bin/sh
# A script to generate the png and avi files

echo "Running $0" "$@"
IP=$1

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
#echo "UUID = $UUID -- args = $args"

DIR="/paasage/$UUID"
cd $DIR/Molecular_Docking

# NOTE: debugging
#echo
#echo $DIR/Molecular_Docking
#ls -lah $DIR/Molecular_Docking
#echo

sh make-image.sh

# If IP of a worker node is given, run this remotely
#ssh -o StrictHostKeyChecking=no -o NumberOfPasswordPrompts=0 \
#        -o LogLevel=quiet -o ConnectTimeout=1 paasage@$IP \
#        "cd $DIR/Molecular_Docking ; sh make-image.sh"

tmp=`ls | grep result_`
mv $tmp $UUID
tar cf $UUID.tar $UUID/
gzip *.tar
mv *.gz /paasage

cd /paasage
rm -rf $DIR/

