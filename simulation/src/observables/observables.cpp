#include"observables.h"

int magnetization(const IsingConfiguration& inConfiguration){
    int magnetization{ 0 };
    const signed char* spins = inConfiguration.get_spins_data();
    for (unsigned int site = 0; site < inConfiguration.get_num_spins(); site++) {magnetization += spins[site]; }
    return magnetization;
}

double energy(const Lattice & inLattice, const J1J2 & inBonds, const IsingConfiguration & inConfiguration){
    double energy{ 0.0 }, delta;
    unsigned int neighbor;
    const unsigned int numSites = inLattice.get_num_sites();
    const unsigned short numNearestNeighbors = inLattice.get_num_nearest_neighbors();
    const unsigned int * allSiteNeighbors = inLattice.get_all_site_neighbors_data();
    const unsigned int * allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();
    const signed char * spins = inConfiguration.get_spins_data();
    const unsigned char * bonds = inBonds.get_bonds_data();
    const double j2 = inBonds.get_j2();
    for(unsigned int site = 0; site < numSites; site++){
        neighbor = site * numNearestNeighbors;
        if(bonds[allEdgeNeighbors[neighbor]]){delta = spins[allSiteNeighbors[neighbor]];}
        else{delta = j2*spins[allSiteNeighbors[neighbor]];}
        neighbor += 2;
        for(unsigned short n = 2; n < numNearestNeighbors; n += 2) {
            if(bonds[allEdgeNeighbors[neighbor]]){delta += spins[allSiteNeighbors[neighbor]];}
            else{delta += j2*spins[allSiteNeighbors[neighbor]];}
            neighbor += 2;
        }
        energy += delta * spins[site];
    }
    return -1. * energy;
}
