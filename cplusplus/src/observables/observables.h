#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "../lattice/lattice.h"
#include "../ising_configuration/ising_configuration.h"
#include "../random_numbers/random_numbers.h"
#include "../interaction/j1j2.h"



int magnetization(const IsingConfiguration & inConfiguration);
double energy(const Lattice & inLattice, const J1J2 & inBonds, const IsingConfiguration & inConfiguration);



#endif