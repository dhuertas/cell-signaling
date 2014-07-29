#!/bin/bash

INI_FILE="./networks/tests/hard-spheres-test-script.ini"
OUT_FILE="temperature-density-fraction.out"
CMD="./cell-signaling -u Cmdenv -f omnetpp.ini -c temperature-density-fraction"

M_PI="3.14159265358979323846"
k="1.3806488*10^-23" #Boltzmann constant
N="1000"
R="0.275/2.0"
m="2.99150736132928*10^-23" # grams

> $OUT_FILE

temp=("20" "40" "60" "80" "100" "120" "140" "160" "180" "200" "220" "240" "260" "280" "300" "320" "340")
defrac=("5" "10" "15" "20" "25" "30" "35" "40" "45" "50" "55" "60")

for T in "${temp[@]}"
do
	
	# Density fraction [0, 55] %
	for DF in "${defrac[@]}"
	do

		echo "[Config temperature-density-fraction]" > $INI_FILE

		echo "HardSpheresTest.molecule[*].vx = normal(0,($k*$T)/($m))" >> $INI_FILE
		echo "HardSpheresTest.molecule[*].vy = normal(0,($k*$T)/($m))" >> $INI_FILE
		echo "HardSpheresTest.molecule[*].vz = normal(0,($k*$T)/($m))" >> $INI_FILE

		echo "HardSpheresTest.spaceSizeX = ($N*(4/3.0)*$M_PI*($R)^3/($DF/100))^(1/3.0)" >> $INI_FILE
		echo "HardSpheresTest.spaceSizeY = ($N*(4/3.0)*$M_PI*($R)^3/($DF/100))^(1/3.0)" >> $INI_FILE
		echo "HardSpheresTest.spaceSizeZ = ($N*(4/3.0)*$M_PI*($R)^3/($DF/100))^(1/3.0)" >> $INI_FILE

		echo "output-vector-file = /home/dani/Workspace/cell-signaling/results/${configname}-${runnumber}-temp-$T-density-fraction-$DF.vec" >> $INI_FILE
		echo "output-scalar-file = /home/dani/Workspace/cell-signaling/results/${configname}-${runnumber}-temp-$T-density-fraction-$DF.sca" >> $INI_FILE

		$CMD >> $OUT_FILE
	done
done