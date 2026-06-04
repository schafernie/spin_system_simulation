#!/bin/bash

seed=0
for T in $(seq 0.70 0.05 1.20); do
    echo "Running T=$T, DisorderSeed=$seed ..."
    ./mc_sim -O 1 -L 32 -D 2 -T $T -J 0.1 -P 0.5 \
        --DisorderSeed $seed --DisorderSample 0 --NrDisorderSamples 10 \
        --EquilibrationTime 2500 --ThermalSeed 0 --ThermalSamples 1000 --TimeBetweenSamples 25
    ((seed++))
done


