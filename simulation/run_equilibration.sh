#!/bin/bash

set -e

NUM_DISORDER_SAMPLES=10

DISORDER_SEED=0
TEMPERATURE=0.908
EQUILIBRATION_TIME=32500


for (( i=0; i < $NUM_DISORDER_SAMPLES ; i++ )); do
    #./simulation/mc_sim -O 0 -L 32 -D 2 -T $TEMPERATURE -J 0.1 -P 0.5 --DisorderSeed $i --EquilibrationTime 32500 
    #echo " $ROOT_DIR/simulation/mc_sim -O 0 -L 32 -D 2 -J 0.1 -P 0.5 -T $TEMPERATURE --DisorderSeed $i --EquilibrationTime $EQUILIBRATION_TIME" 
    cmd="./mc_sim -O 0 -L 32 -D 2 -J 0.1 -P 0.5 -T $TEMPERATURE \
    --DisorderSeed $DISORDER_SEED --DisorderSample $i --EquilibrationTime $EQUILIBRATION_TIME" 
    echo $cmd
    $cmd
done
