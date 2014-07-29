#!/bin/bash

me=`basename $0`
cmd="./cell-signaling -u Cmdenv"
omnet_ini="./omnetpp.ini"

M_PI="3.14159265358979323846"

rt_vs_number_of_particles() {

  ini_file="./networks/tests/rt-vs-number-of-particles.ini"
  out_file="./results/rt-vs-number-of-particles.out"
  config="rt-vs-number-of-particles"

  > $out_file

  echo "### $config ###" >> $out_file

  vf="36"

  for sl in {2..8..6}
  do
    for mode in {1..2}
    do 
      for i in {25000..60000..5000}
      do
        > $ini_file
        echo "[Config $config]" >> $ini_file
 
        echo "HardSpheresTest.numberOfInitialMolecules = $i" >> $ini_file
        echo "HardSpheresTest.spaceSizeX = ($i*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file
        echo "HardSpheresTest.spaceSizeY = ($i*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file
        echo "HardSpheresTest.spaceSizeZ = ($i*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file

        echo "HardSpheresTest.manager.mode = $mode" >> $ini_file
        echo "HardSpheresTest.manager.spaceCellSize = $sl" >> $ini_file

        echo "HardSpheresTest.molecule[*].listRadius = 2" >> $ini_file
        echo "HardSpheresTest.molecule[*].refreshListRadius = 2" >> $ini_file
		
        #ts=`date +%Y%m%d_%H%M%S`
        #echo "output-vector-file = \${resultdir}/$ts-\${configname}-\${runnumber}-$i.vec" >> $ini_file
        #echo "output-scalar-file = \${resultdir}/$ts-\${configname}-\${runnumber}-$i.sca" >> $ini_file

        echo "Side length $sl, mode $mode, number of particles $i" >> $out_file
        #$cmd -c $config -f $omnet_ini | grep -i "100%\|Error" >> $out_file
        $cmd -c $config -f $omnet_ini >> $out_file
      done
    done
  done
}

rt_vs_volume_fraction() {

  ini_file="./networks/tests/rt-vs-volume-fraction.ini"
  out_file="./results/rt-vs-volume-fraction.out"
  config="rt-vs-volume-fraction"

  > $out_file

  echo "### $config ###" >> $out_file
 
  for vf in {2..40..2}
  do
    > $ini_file
    echo "[Config $config]" >> $ini_file
    echo "HardSpheresTest.numberOfInitialMolecules = 10000" >> $ini_file

    echo "HardSpheresTest.spaceSizeX = (10000*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file
    echo "HardSpheresTest.spaceSizeY = (10000*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file
    echo "HardSpheresTest.spaceSizeZ = (10000*(4/3.0)*$M_PI/($vf/100))^(1/3.0)" >> $ini_file
    
    echo "Volume fraction $vf" >> $out_file
    $cmd -c $config -f $omnet_ini | grep -i "100%\|Error" >> $out_file
  done
}

rt_vs_space_cell_length_vs_verlet_radius() {

  ini_file="./networks/tests/rt-vs-space-cell-length-vs-verlet-radius.ini"
  out_file="./results/rt-vs-space-cell-length-vs-verlet-radius.out"
  config="rt-vs-space-cell-length-vs-verlet-radius"

  #> $out_file

  #echo "### $config ###" >> $out_file

  VR=( "1.1" "1.2" "1.3" "1.4" "1.5" "1.6" "1.7" "1.8" "1.9" )
 
  for sl in {4..8..2}
  do
    for vr in "${VR[@]}"
    do
      > $ini_file
      echo "[Config $config]" >> $ini_file
      echo "HardSpheresTest.numberOfInitialMolecules = 10000" >> $ini_file

      echo "HardSpheresTest.spaceSizeX = (10000*(4/3.0)*$M_PI/(38/100))^(1/3.0)" >> $ini_file
      echo "HardSpheresTest.spaceSizeY = (10000*(4/3.0)*$M_PI/(38/100))^(1/3.0)" >> $ini_file
      echo "HardSpheresTest.spaceSizeZ = (10000*(4/3.0)*$M_PI/(38/100))^(1/3.0)" >> $ini_file

      echo "HardSpheresTest.manager.spaceCellSize = $sl" >> $ini_file
      echo "HardSpheresTest.molecule[*].listRadius = $vr" >> $ini_file
      echo "HardSpheresTest.molecule[*].refreshListRadius = $vr" >> $ini_file
    
      echo "Side length $sl, verlet radius $vr" >> $out_file
      $cmd -c $config -f $omnet_ini | grep -i "100%\|Error" >> $out_file
    done
  done
}

diffusion_pulse() {

  ini_file="./networks/tests/test-emitter-brownian-motion.ini"
  out_file="./results/test-emitter-brownian-motion.out"
  config="test-emitter-brownian-motion"

  > $out_file
  for D in {1..10..1}
  do
    > $ini_file
    echo "[Config $config]" >> $ini_file
    echo "domain.cell[0].emitter.emissionDiffusion = $D" >> $ini_file
    
    ts=`date +%Y%m%d_%H%M%S`
    echo "output-vector-file = \${resultdir}/$ts-\${configname}-\${runnumber}-diffusion-$D.vec" >> $ini_file
    echo "output-scalar-file = \${resultdir}/$ts-\${configname}-\${runnumber}-diffusion-$D.sca" >> $ini_file

    echo "Diffusion $D" >> $out_file
    $cmd -c $config -f $omnet_ini >> $out_file
  done

}

diffusion_pulse_with_seed() {

  results="/home/dani/Workspace/cell-signaling/results"
  ini_file="./networks/tests/test-emitter-brownian-motion.ini"
  out_file="./results/test-emitter-brownian-motion.out"
  config="test-emitter-brownian-motion"

  > $out_file
  for i in {1..30}
  do
    > $ini_file
    seed=`od -vAn -N4 -tu4 < /dev/urandom`
    echo "[Config $config]" >> $ini_file
    echo "seed-0-mt = $seed" >> $ini_file
    echo "domain.cell[0].emitter.emissionDiffusion = 1" >> $ini_file

    ts=`date +%Y%m%d_%H%M%S`
    echo "output-vector-file = $results/$ts-\${configname}-\${runnumber}-seed-$seed.vec" >> $ini_file
    echo "output-scalar-file = $results/$ts-\${configname}-\${runnumber}-seed-$seed.sca" >> $ini_file

    echo "Seed $seed" >> $out_file
    $cmd -c $config -f $omnet_ini >> $out_file
  done
}

#rt_vs_number_of_particles

#rt_vs_volume_fraction

#rt_vs_space_cell_length_vs_verlet_radius

#diffusion_pulse

diffusion_pulse_with_seed
