#!/bin/sh
# Copyright 2013 University of Stuttgart, Germany
# Author: Anthony Sulistio (HLRS)
#
# A script to set up NFS and create hostfile.txt that lists IPs of worker nodes.
# This script needs to be executed from the MPI master node
# Usage:   ./setup_NFS.sh NFS_IP [start_IP] [end_IP]
#
# Example: ./setup_NFS.sh 109.231.122.193  170  220 -> will create an NFS server 
# on 109.231.122.193 and scan IP from 109.231.122.170 to 109.231.122.220 for
# list of MPI worker nodes. 
# The mounted NFS directory will be in 109.231.122.193:/paasage
#
# Note: for this script to work properly, here are the steps:
# - generate ssh key : ssh-keygen -t rsa
# - put id_rsa.pub into /root/.ssh/authorized_keys and /home/paasage/.ssh/authorized_keys
# - put id_rsa into /root/.ssh/id_rsa and /home/paasage/.ssh/id_rsa

## running the hostfile on MPI: -hostfile hostfile.txt
#mpirun -np $CPU -hostfile hostfile.txt ./main --cpu $CPU -N $NUM  ...

## mount NFS
#mount -t nfs $STORAGE_IP:/paasage /paasage

IPADDR=`ifconfig eth0 | grep -i bcast | gawk -F" " '{ print $2 }' | gawk -F":" '{ print $2 }'`

# set up the NFS server
if [ -n "$1" ]; then
    NFS=$1
    NODEIP=`echo $NFS | sed 's/\./ /g' | awk {'print $1"."$2"."$3'}`
    CMD="echo \"/paasage $NODEIP.0/24(rw,sync,no_root_squash,no_subtree_check)\" >> /etc/exports"
    echo $CMD
    ssh -o StrictHostKeyChecking=no -o NumberOfPasswordPrompts=0 \
        -o LogLevel=quiet -o ConnectTimeout=1 root@$NFS \
        "$CMD ; /etc/init.d/nfs-kernel-server restart ; /etc/init.d/lighttpd stop"
fi

NODEIP=`echo $IPADDR | sed 's/\./ /g' | awk {'print $1"."$2"."$3'}`
TXT=hostfile.txt

## starting IP address
index=1
if [ -n "$2" ]; then
    index=$2
fi

## ending IP address
endIP=254
if [ -n "$3" ]; then
    endIP=$3
fi

echo "Server IP:" $IPADDR
echo "NFS IP:" $NFS
echo "Searching worker nodes from:" $NODEIP.$index "to"  $NODEIP.$endIP
echo

# remove the files first
rm -f $TXT
rm -f ~/.ssh/known_hosts

while [ $index -le $endIP ]; do
    ADDR=`echo $NODEIP"."$index`
    echo "Scanning $ADDR ..."
    index=`expr $index + 1`

    if [ "$ADDR" = "$NFS" ]; then
        echo "-- $NFS is the NFS server. Ignore."
        continue
    fi

    if [ "$ADDR" = "$IPADDR" ]; then
        echo "-- $IPADDR is the MPI master node. Ignore."
        mount -t nfs $NFS:/paasage /paasage
        continue
    fi

    # connect the worker node and mount the NFS directory
    NODEPROCS=`ssh -o StrictHostKeyChecking=no -o NumberOfPasswordPrompts=0 \
        -o LogLevel=quiet -o ConnectTimeout=1 root@$ADDR \
        "mount -t nfs $NFS:/paasage /paasage ; /etc/init.d/lighttpd stop > /dev/null ; cat /proc/cpuinfo | grep processor | wc -l"`

    if [ -n "$NODEPROCS" ]; then
        echo "-- $ADDR has $NODEPROCS cores"
        echo "$ADDR  slots=$NODEPROCS" >> $TXT 
    fi
done

chmod 666 $TXT
chown paasage:paasage $TXT

echo
echo "Check $TXT for completeness and correctness."
echo "Then, copy $TXT into the directory where a program will be executed"
echo

