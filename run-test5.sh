#!/bin/bash

ME=`basename $0`

CMD="./CellSignaling -u Cmdenv"

OMNET_INI="./omnetpp.ini"

NETWORK="test5"
NETWORKS_FOLDER="./networks"

NETWORK_INI="$NETWORKS_FOLDER/$NETWORK.ini"

OUTPUT_FOLDER="./results"
OUTPUT_FILE="$OUTPUT_FOLDER/$ME.out"

# Remove previous file
rm -rf $OUTPUT_FILE

M_PI="3.14159265358979323846"
VD="0.00"
VD_STEP="0.01"

# $1 The number of initial particles per axis
# $2 The volume density
write_ini() {
	echo "[Config $NETWORK]" > $NETWORK_INI
	echo "$NETWORK.numberOfInitialMolecules = \${N=$1}^3" >> $NETWORK_INI
	echo "$NETWORK.spaceSizeX = (\${N}^3*(4/3.0)*$M_PI/$2)^(1/3.0)" >> $NETWORK_INI
	echo "$NETWORK.spaceSizeY = (\${N}^3*(4/3.0)*$M_PI/$2)^(1/3.0)" >> $NETWORK_INI
	echo "$NETWORK.spaceSizeZ = (\${N}^3*(4/3.0)*$M_PI/$2)^(1/3.0)" >> $NETWORK_INI
}

VD="0.01"
for i in {10..30..2}
do
	write_ini $i $VD
	$CMD -u Cmdenv -c $NETWORK -f $OMNET_INI | grep "100%" | awk '{print $6}' >> $OUTPUT_FILE
done

# for i in {1..50}
# do
#	VD=`echo "$VD $VD_STEP" | awk '{printf "%f", $1 + $2}'`
#	write_ini 10 $VD
#	$CMD -u Cmdenv -c test4 -f $OMNET_INI | grep "100%" | awk '{print $6}' >> $OUTPUT_FILE
# done
