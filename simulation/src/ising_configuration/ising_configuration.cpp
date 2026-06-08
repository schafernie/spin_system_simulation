#include "ising_configuration.h"


IsingConfiguration::IsingConfiguration(const Lattice& inLattice)
{
    numSpins = inLattice.get_num_sites();
    spins.resize(numSpins);
    set_spins(1);
}


void IsingConfiguration::set_spins(const int inS){
    signed char s0;
    if(inS >= 0){
        s0 = 1;
    }
    else{
        s0 = -1;
    }   
    for(unsigned int c=0; c<numSpins; c++){
        spins[c] = s0;
    }
}


void IsingConfiguration::set_spins_rand(const double inP, RandomNumbers& inRandomNumbers){
    std::vector<unsigned char> bernoullis = inRandomNumbers.bernoulli_vector(inP, numSpins); 
    for (unsigned int c = 0; c < numSpins; c++) {
        if (bernoullis[c]) {
            spins[c] = 1;
        }
        else {
            spins[c] = -1;
        }
    }
}


void IsingConfiguration::print(const Lattice& InLattice)const{
    const unsigned short dimension = InLattice.get_dimension();
    const unsigned int linearSize = InLattice.get_linear_size();
    const unsigned int * dimensionFactors = InLattice.get_dimension_factors_data();  
    
    std::cout << "\nIsing configuration:";
    std::cout << "\ndimension: " << dimension;
    std::cout << "\nlinear size: " << linearSize;
    std::cout << "\n";
    std::vector<unsigned int> LatticeIndex(dimension, 0);
    unsigned int site;

    std::cout << std::left;
    for (unsigned int c0 = 0; c0 < numSpins; c0++) {
        site = 0;
        for (unsigned short c1 = 0; c1 < dimension; c1++) {
            site += dimensionFactors[c1] * LatticeIndex[c1];
        }
        std::cout.width(4);
        std::cout << static_cast<int>(spins[site]);
        LatticeIndex[0]++;
        for (unsigned short c2 = 1; c2 < dimension; c2++) {
            if ((LatticeIndex[c2 - 1] > 0) && ((LatticeIndex[c2 - 1] % linearSize) == 0)) {
                LatticeIndex[c2 - 1] = 0;
                LatticeIndex[c2]++;
                std::cout << "\n";
            }
        }
    }
    std::cout << "\n";
}
