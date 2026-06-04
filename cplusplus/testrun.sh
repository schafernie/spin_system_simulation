#equilibration
./mc_sim -O 0 -L 32 -D 2 -T 0.90 -J 0.1 -P 0.5 --DisorderSeed 0 --EquilibrationTime 25000 

#disorder sampling
./mc_sim -O 1 -L 32 -D 2 -T 0.90 -J 0.1 -P 0.5 --DisorderSeed 0 --NrDisorderSamples 10 --EquilibrationTime 1000 --ThermalSeed 0 --ThermalSamples 1000 --TimeBetweenSamples 30 

#sampling single sample 
./mc_sim -O 1 -L 32 -D 2 -T 0.90 -J 0.1 -P 0.5 --DisorderSeed 0 --NrDisorderSamples 10 --EquilibrationTime 1000 --ThermalSeed 0 --ThermalSamples 1000 --TimeBetweenSamples 30 --FixedDisorder 
