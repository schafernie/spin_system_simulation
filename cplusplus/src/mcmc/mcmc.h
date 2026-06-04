#ifndef MCMC_H
#define MCMC_H


#include "../lattice/lattice.h"
#include "../clusters/clusters.h"
#include "../ising_configuration/ising_configuration.h"
#include "../random_numbers/random_numbers.h"
#include "../interaction/j1j2.h"

  
//single spin-flip moves 
void mcmc_single_spin_flip(const Lattice & inLattice, const J1J2 & inBonds, const double inBeta, IsingConfiguration & InConfiguration,
    RandomNumbers & inRandomNumbers, const unsigned int inNumSteps = 1);

//Swendsen-Wang moves 
void mcmc_swendsen_wang(const Lattice & inLattice, const J1J2 & inBonds, const double  inBeta, IsingConfiguration & inConfiguration, 
    Clusters & inClusters, RandomNumbers & inRandomNumbers,  const unsigned int inNumSteps = 1);

#endif 
