#include"mcmc.h"


void mcmc_single_spin_flip(const Lattice & inLattice, const J1J2 & inBonds, const double inBeta, IsingConfiguration & inConfiguration,
    RandomNumbers & inRandomNumbers, const unsigned int inNumSteps){
    const unsigned int numSites = inLattice.get_num_sites();
    const unsigned short numNearestNeighbors = inLattice.get_num_nearest_neighbors();
    const unsigned int* allSiteNeighbors = inLattice.get_all_site_neighbors_data();
    const unsigned int* allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();
    const unsigned char * bonds = inBonds.get_bonds_data();
    const double j2 = inBonds.get_j2();

    signed char * spins = inConfiguration.get_spins_data();
    unsigned int index = inRandomNumbers.get_index();
    double * randomBuffer01 = inRandomNumbers.get_random_buffer01_data();
    int * randomIntegers = inRandomNumbers.get_random_integers_data(); 

    unsigned int site, neighbor;
    unsigned short r;
    double delta;
    for (unsigned int step = 0; step < inNumSteps; step++) {
        inRandomNumbers.fill_random_buffer_01(index, numSites);
        inRandomNumbers.generate_uniform_random_integers(numSites, 0, numSites-1);
        for(unsigned int c = 0; c < numSites; c++) {
            site = randomIntegers[c];
            neighbor = site * numNearestNeighbors;
            if (bonds[allEdgeNeighbors[neighbor]]) { delta = spins[allSiteNeighbors[neighbor]]; }
            else { delta = j2 * spins[allSiteNeighbors[neighbor]]; }
            neighbor++;
            for (r = 1; r < numNearestNeighbors; r++) {
                if (bonds[allEdgeNeighbors[neighbor]]) {
                    delta += spins[allSiteNeighbors[neighbor]];
                }
                else{ delta += j2 * spins[allSiteNeighbors[neighbor]]; }
                neighbor++;
            }
            delta *= spins[site];
            if (delta <= 0) { spins[site] *= -1; }
            else if (randomBuffer01[index++] < std::exp(-2.0 * inBeta * delta)) { spins[site] *= -1; }
        }
    }
    inRandomNumbers.set_index(index);
}


void mcmc_swendsen_wang(const Lattice & inLattice, const J1J2 & inBonds, const double inBeta,
     IsingConfiguration & inConfiguration, Clusters & inClusters, RandomNumbers & inRandomNumbers, const unsigned int inNumSteps)
    {
    const unsigned int * allSiteNeighbors = inLattice.get_all_site_neighbors_data();
    const unsigned int * allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();     
    const unsigned int * allAdjacentSites = inLattice.get_all_adjacent_sites_data();
    const unsigned char * bonds = inBonds.get_bonds_data();

    const unsigned int numSites = inLattice.get_num_sites();
    const unsigned int numEdges = inLattice.get_num_edges();
    const unsigned short numNearestNeighbors = inLattice.get_num_nearest_neighbors();

    signed char * spins = inConfiguration.get_spins_data();
    unsigned int index = inRandomNumbers.get_index();
    double * randomBuffer01 = inRandomNumbers.get_random_buffer01_data();

    int * clusterIds = inClusters.get_cluster_ids_data();
    unsigned char * occupiedEdges = inClusters.get_occupied_edges_data();
    int * stack = inClusters.get_stack_data();

    const double occProb1 = 1.0 - std::exp(-2.0 * inBeta);
    const double occProb2 = 1.0 - std::exp(-2.0 * std::abs(inBonds.get_j2()) * inBeta);

    unsigned int stackCounter;
    unsigned int neighborEdge, neighborSite, currentSite, neighbor, currentEdge, site, currentClusterId;
    bool doFlip;

    for (unsigned int n = 0; n < inNumSteps; n++) {
        inRandomNumbers.fill_random_buffer_01(index, 2 * numEdges);

        for (site = 0; site < numSites; site++) clusterIds[site] = -1;

        for (currentEdge = 0; currentEdge < numEdges; currentEdge++) {
            if ((spins[allAdjacentSites[2 * currentEdge]] * spins[allAdjacentSites[2 * currentEdge + 1]]) > 0) {
                const double occProb = bonds[currentEdge] ? occProb1 : occProb2;
                occupiedEdges[currentEdge] = (randomBuffer01[index++] < occProb) ? 1 : 0;
            }
            else { occupiedEdges[currentEdge] = 0; }
        }
        currentClusterId = 0;
        for (site = 0; site < numSites; site++) {
            if (clusterIds[site] < 0) {
                if(randomBuffer01[index++] < 0.5){doFlip = true; }
                else {doFlip = false; }
                clusterIds[site] = currentClusterId;
                stack[0] = site;
                stackCounter = 1;
                while (stackCounter > 0) {
                    currentSite = stack[--stackCounter];
                    if (doFlip) { spins[currentSite] *= -1; }
                    neighbor = currentSite * numNearestNeighbors;
                    for (unsigned int r = neighbor; r < (neighbor + numNearestNeighbors); r++) {
                        neighborSite = allSiteNeighbors[r];
                        neighborEdge = allEdgeNeighbors[r];
                        if (occupiedEdges[neighborEdge] && (clusterIds[neighborSite] < 0)) {
                            clusterIds[neighborSite] = currentClusterId;
                            stack[stackCounter++] = neighborSite;
                        }
                    }
                }
               currentClusterId++;
            }
        }
    }
    inRandomNumbers.set_index(index);
}