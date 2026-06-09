#!/bin/bash

set -e

NUM_DISORDER_SAMPLES=10

DISORDER_SEED=0
TEMPERATURE=0.905883
EQUILIBRATION_TIME=32500


for (( i=0; i < $NUM_DISORDER_SAMPLES ; i++ )); do
    cmd="./mc_sim -O 0 -L 32 -D 2 -J 0.1 -P 0.5 -T $TEMPERATURE \
    --DisorderSeed $DISORDER_SEED --DisorderSample $i --EquilibrationTime $EQUILIBRATION_TIME" 
    echo $cmd
    $cmd
done
