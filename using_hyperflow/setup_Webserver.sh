#!/bin/sh
# Copyright 2013 University of Stuttgart, Germany
# Author: Anthony Sulistio (HLRS)
#
# A script to set up a webserver
# Usage: ./setup_Webserver.sh

#apt-get install lighttpd

FILE="/etc/lighttpd/lighttpd.conf"
count=`grep -i paasage $FILE | wc -l`
if [ $count -eq 0 ]; then
    #echo "Updating $FILE"
    echo 'alias.url += ( "/paasage" => "/paasage/" )' >> $FILE
    echo '$HTTP["url"] =~ "^/paasage($|/)" { server.dir-listing = "enable" }' >> $FILE
fi

/etc/init.d/lighttpd restart

echo "Installing the webserver is done."
echo "The URL for PaaSage directory is http://127.0.0.1/paasage"
echo

