#!/bin/bash

ME=`basename $0`
CMD="./CellSignaling -u Cmdenv"
OMNET_INI="./omnetpp.ini"

OUTPUT_FOLDER="./results"
OUTPUT_FILE="$OUTPUT_FOLDER/$ME.out"

NETWORKS_FOLDER="./networks"

M_PI="3.14159265358979323846"

# --------------------------------------------------------------------------- #
rt_vs_number_of_particles() {

	echo "#############################" >> $OUTPUT_FILE
	echo "# rt-vs-number-of-particles #" >> $OUTPUT_FILE
	echo "#############################" >> $OUTPUT_FILE

	CONFIG="rt-vs-number-of-particles"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-number-of-particles.ini"

	for i in {10000..100000..10000}
	do
		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.numberOfInitialMolecules = $i" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.sca" >> $INI_FILE

		echo "Number of particles $i" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
rt_vs_particle_distribution() {
	echo "###############################" >> $OUTPUT_FILE
	echo "# rt-vs-particle-distribution #" >> $OUTPUT_FILE
	echo "###############################" >> $OUTPUT_FILE

	CONFIG="rt-vs-particle-distribution"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-particle-distribution.ini"

	DISTRIB=( "uniform" "sphere" "cube" "highdensity")

	for i in ${DISTRIB[@]}
	do
		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.manager.particleDistribution = \"$i\"" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.sca" >> $INI_FILE

		echo "Distribution $i" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
rt_vs_space_cell_length() {
	echo "###########################" >> $OUTPUT_FILE
	echo "# rt-vs-space-cell-length #" >> $OUTPUT_FILE
	echo "###########################" >> $OUTPUT_FILE

	CONFIG="rt-vs-space-cell-length"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-space-cell-length.ini"

	VD="0.15"

	for i in {2..10..1}
	do
		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.manager.spaceCellSize = $i" >> $INI_FILE

		echo "hsa.numberOfInitialMolecules = \${N=10}^3" >> $INI_FILE
		echo "hsa.spaceSizeX = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeY = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeZ = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.sca" >> $INI_FILE

		echo "Space cell size $i" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
rt_vs_verlet_list_radius() {
	echo "############################" >> $OUTPUT_FILE
	echo "# rt-vs-verlet-list-radius #" >> $OUTPUT_FILE
	echo "############################" >> $OUTPUT_FILE

	CONFIG="rt-vs-verlet-list-radius"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-verlet-list-radius.ini"

	VD="0.15"

	for i in {12..30..2}
	do
		LR=`echo "$i 10" | awk '{printf "%f", $1/$2}'`

		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.molecules[*].listRadius = $LR" >> $INI_FILE

		echo "hsa.numberOfInitialMolecules = \${N=10}^3" >> $INI_FILE
		echo "hsa.spaceSizeX = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeY = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeZ = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}.sca" >> $INI_FILE

		echo "List radius $LR" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
rt_vs_volume_density() {
	echo "########################" >> $OUTPUT_FILE
	echo "# rt-vs-volume-density #" >> $OUTPUT_FILE
	echo "########################" >> $OUTPUT_FILE

	CONFIG="rt-vs-volume-density"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-volume-density.ini"

	VD="0.00"
	VD_STEP="0.05"

	for i in {1..10}
	do
		VD=`echo "$VD $VD_STEP" | awk '{printf "%f", $1 + $2}'`

		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.numberOfInitialMolecules = \${N=10}^3" >> $INI_FILE
		echo "hsa.spaceSizeX = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeY = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
		echo "hsa.spaceSizeZ = (\${N}^3*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$VD.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$VD.sca" >> $INI_FILE

		echo "Volume density $VD" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
number_collisions_over_simtime() {
	echo "##################################" >> $OUTPUT_FILE
	echo "# number-collisions-over-simtime #" >> $OUTPUT_FILE
	echo "##################################" >> $OUTPUT_FILE

	CONFIG="number-collisions-over-simtime"
	INI_FILE="$NETWORKS_FOLDER/number-collisions-over-simtime.ini"

	DISTRIB=( "uniform" "sphere" "cube" "highdensity" )
	for i in ${DISTRIB[@]}
	do
		echo "[Config $CONFIG]" > $INI_FILE
		echo "hsa.manager.particleDistribution 	= \"$i\"" >> $INI_FILE

		RUN_TS=`date +%Y%m%d_%H%M%S`
		echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.vec" >> $INI_FILE
		echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$i.sca" >> $INI_FILE

		echo "Distribution $i" >> $OUTPUT_FILE
		$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
	done
}

# --------------------------------------------------------------------------- #
space_position_over_simtime() {
	echo "###############################" >> $OUTPUT_FILE
	echo "# space-position-over-simtime #" >> $OUTPUT_FILE
	echo "###############################" >> $OUTPUT_FILE

	CONFIG="space-position-over-simtime"
	INI_FILE="$NETWORKS_FOLDER/space-position-over-simtime.ini"

	echo "[Config $CONFIG]" > $INI_FILE

	RUN_TS=`date +%Y%m%d_%H%M%S`
	echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}.vec" >> $INI_FILE
	echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}.sca" >> $INI_FILE

	$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE

}

# --------------------------------------------------------------------------- #
rt_vs_space_cell_length_vs_verlet_radius() {
	echo "############################################" >> $OUTPUT_FILE
	echo "# rt-vs-space-cell-length-vs-verlet-radius #" >> $OUTPUT_FILE
	echo "############################################" >> $OUTPUT_FILE

	CONFIG="rt-vs-space-cell-length-vs-verlet-radius"
	INI_FILE="$NETWORKS_FOLDER/rt-vs-space-cell-length-vs-verlet-radius.ini"

	VDS=( "0.15" "0.35")

	for VD in ${VDS[@]} # Volume density
	do
		for CS in {2..10} # Space cell size (particle radius = 1, diameter = 2)
		do
			VLR=( "1.2" "1.4" "1.6" "1.8" "2.0" "2.2" "2.4" "2.6" "2.8" "3.0" ) # Verlet list radius
			for LR in ${VLR[@]}
			do
				echo "[Config $CONFIG]" > $INI_FILE
				echo "hsa.manager.spaceCellSize = $CS" >> $INI_FILE
				echo "hsa.numberOfInitialMolecules = \${N=10000}" >> $INI_FILE
				echo "hsa.spaceSizeX = (\${N}*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
				echo "hsa.spaceSizeY = (\${N}*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
				echo "hsa.spaceSizeZ = (\${N}*(4/3.0)*$M_PI/$VD)^(1/3.0)" >> $INI_FILE
				echo "hsa.molecules[*].listRadius = $LR" >> $INI_FILE

				RUN_TS=`date +%Y%m%d_%H%M%S`
				echo "output-vector-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$VD-$CS-$LR.vec" >> $INI_FILE
				echo "output-scalar-file = ../\${resultdir}/$RUN_TS-\${configname}-\${runnumber}-$VD-$CS-$LR.sca" >> $INI_FILE

				echo "Volume density $VD, space cell size $CS, Verlet list radius $LR" >> $OUTPUT_FILE
				$CMD -c $CONFIG -f $OMNET_INI | grep -i "100%\|Error" >> $OUTPUT_FILE
			done
		done
	done
}

rm $OUTPUT_FILE

# rt_vs_number_of_particles
# rt_vs_particle_distribution
# rt_vs_space_cell_length
# rt_vs_verlet_list_radius
# rt_vs_volume_density
# number_collisions_over_simtime
# space_position_over_simtime
rt_vs_space_cell_length_vs_verlet_radius