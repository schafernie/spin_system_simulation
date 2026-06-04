#ifndef PROTOCOL_H
#define PROTOCOL_H

#include <iostream>


#include"../lattice/lattice.h"
#include"../interaction/j1j2.h"
#include"../ising_configuration/ising_configuration.h"
#include"../random_numbers/random_numbers.h"
#include"../observables/observables.h"
#include"../mcmc/mcmc.h"
#include"../statistics/statistics.h"

using namespace std;

struct Parameters {
    unsigned short dimension{2};
    unsigned int linearSize{4};
    double temperature{1.};
    double p2{0.5};
    double j2{0.1};
    int disorderSeed{0};
    int disorderSample{0};
    int thermalSeed{0};
    unsigned int thermalEquilibrationTime{1000};
    unsigned int thermalSamples{100};
    unsigned int thermalTimeBetweenSamples{10};
    unsigned int nrDisorderSamples{1};
    bool fixedDisorder{false};
};



void j1j2_equilibration(const Parameters & );
void j1j2_disorder_sampling(const Parameters &);
void debugging(const Parameters & );



#endif

