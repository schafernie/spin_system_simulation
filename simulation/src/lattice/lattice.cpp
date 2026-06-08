#include"lattice.h"


Lattice::Lattice(){}


Lattice::Lattice(const unsigned short inDimension, const unsigned int inLinearSize)
{
    dimension          = inDimension;
    linearSize         = inLinearSize;
    numNearestNeighbors  = 2 * dimension;
    numSites           = static_cast<unsigned int>(std::round(std::pow(linearSize, dimension)));
    numEdges           = dimension * numSites;
    numAllNeighbors    = numNearestNeighbors * numSites;
    numAllAdjacentSites = 2 * numEdges;
    numAllCoordinates  = dimension * numSites;

    allSiteNeighbors.resize(numAllNeighbors);
    allEdgeNeighbors.resize(numAllNeighbors);
    allAdjacentSites.resize(numAllAdjacentSites);
    allCoordinates.resize(numAllCoordinates);
    dimensionFactors.resize(dimension);

    std::vector<unsigned int> latticeIndex(dimension, 0);
    std::vector<unsigned int> currentLatticeIndex(dimension);

    // Compute dimension factors (powers of linearSize)
    for (unsigned short dimIdx = 0; dimIdx < dimension; dimIdx++) {
        unsigned int factor = 1;
        for (unsigned int d = 0; d < dimIdx; d++) {
            factor *= linearSize;
        }
        dimensionFactors[dimIdx] = factor;
    }

    // Compute site neighbors (periodic boundary conditions)
    unsigned int neighborIdx = 0;
    while (neighborIdx < numAllNeighbors) {
        // Compute site neighbors
        for (unsigned short dimIdx = 0; dimIdx < dimension; dimIdx++) {
            // Copy current lattice index
            for (unsigned short d = 0; d < dimension; d++) {
                currentLatticeIndex[d] = latticeIndex[d];
            }

            // Neighbor in plus direction
            currentLatticeIndex[dimIdx] = latticeIndex[dimIdx] + 1;
            if (currentLatticeIndex[dimIdx] > linearSize - 1) {
                currentLatticeIndex[dimIdx] -= linearSize;
            }
            allSiteNeighbors[neighborIdx] = 0;
            for (unsigned short d = 0; d < dimension; d++) {
                allSiteNeighbors[neighborIdx] += dimensionFactors[d] * currentLatticeIndex[d];
            }
            neighborIdx++;

            // Neighbor in minus direction
            currentLatticeIndex[dimIdx] =   
            (latticeIndex[dimIdx] == 0) ? linearSize - 1 : latticeIndex[dimIdx] - 1;
            allSiteNeighbors[neighborIdx] = 0;
            for (unsigned short d = 0; d < dimension; d++) {
                allSiteNeighbors[neighborIdx] += dimensionFactors[d] * currentLatticeIndex[d];
            }
            neighborIdx++;
        }

        // Increment lattice index
        latticeIndex[0]++;
        for (unsigned short dimIdx = 1; dimIdx < dimension; dimIdx++) {
            if ((latticeIndex[dimIdx - 1] > 0) && ((latticeIndex[dimIdx - 1] % linearSize) == 0)) {
                latticeIndex[dimIdx - 1] = 0;
                latticeIndex[dimIdx]++;
            }
        }
    }

    // Compute edge neighbors
    unsigned int edgeIdx = 0;
    unsigned int siteIdx = 0;
    while (edgeIdx < numAllNeighbors) {
        for (unsigned int dimIdx = 0; dimIdx < dimension; dimIdx++) {
            allEdgeNeighbors[edgeIdx] = dimension * siteIdx + dimIdx;           // plus direction
            edgeIdx++;
            allEdgeNeighbors[edgeIdx] = dimension * allSiteNeighbors[numNearestNeighbors * siteIdx + 2 * dimIdx + 1] + dimIdx; // minus direction
            edgeIdx++;
        }
        siteIdx++;
    }

    // Compute adjacent sites
    unsigned int adjIdx = 0;
    for (unsigned int siteIdx = 0; siteIdx < numSites; siteIdx++) {
        for (unsigned short dimIdx = 0; dimIdx < dimension; dimIdx++) {
            allAdjacentSites[adjIdx++] = siteIdx;
            allAdjacentSites[adjIdx++] = allSiteNeighbors[2 * dimension * siteIdx + 2 * dimIdx];
        }
    }

    // Compute coordinates
    std::fill(latticeIndex.begin(), latticeIndex.end(), 0);
    for (unsigned int siteIdx = 0; siteIdx < numSites; siteIdx++) {
        for (unsigned short dimIdx = 0; dimIdx < dimension; dimIdx++) {
            allCoordinates[dimension * siteIdx + dimIdx] = latticeIndex[dimIdx];
        }

        // Increment lattice index
        latticeIndex[0]++;
        for (unsigned short dimIdx = 1; dimIdx < dimension; dimIdx++) {
            if ((latticeIndex[dimIdx - 1] > 0) && ((latticeIndex[dimIdx - 1] % linearSize) == 0)) {
                latticeIndex[dimIdx - 1] = 0;
                latticeIndex[dimIdx]++;
            }
        }
    }
}


Lattice::~Lattice(){
}



void Lattice::print()const{
    std::cout << "\n****periodic symmetric lattice****";
    std::cout << "\ndimension: " << dimension; 
    std::cout << "\nlinear size: " << linearSize;
    std::cout << "\ndimension factors: ";
    for(unsigned short c=0; c<dimension; c++){
        std::cout<<dimensionFactors[c] <<", ";
    }
    std::cout << "\n****sites****\n";
    std::vector<unsigned int> latticeIndex(dimension, 0);
    unsigned int site;
    std::cout << std::left;
    for (unsigned int c0 = 0; c0 < numSites; c0++) {
        site = 0;
        for (unsigned short c1 = 0; c1 < dimension; c1++) {
            site += dimensionFactors[c1] * latticeIndex[c1];
        }
        std::cout.width(4);
        std::cout << site;
        latticeIndex[0]++;
        for (unsigned short c2 = 1; c2 < dimension; c2++) {
            if ((latticeIndex[c2 - 1] > 0) && ((latticeIndex[c2 - 1] % linearSize) == 0)) {
                latticeIndex[c2 - 1] = 0;
                latticeIndex[c2]++;
                std::cout << "\n";
            }
        }
    }
    std::cout<<"\n";
    std::fill(latticeIndex.begin(), latticeIndex.end(),0);
    for (unsigned int c0 = 0; c0 < numSites; c0++) {
        site = 0;
        for (unsigned short c1 = 0; c1 < dimension; c1++) {
            site += dimensionFactors[c1] * latticeIndex[c1];
        }
        std::cout << "[ ";
        for(short c1 = static_cast<short>(dimension - 1); c1 >= 0; c1--) {std::cout << latticeIndex[c1] << " ";}
        std::cout << "] = ";
        std::cout.width(4);
        std::cout << site;
        latticeIndex[0]++;
        for (unsigned short c2 = 1; c2 < dimension; c2++) {
            if ((latticeIndex[c2 - 1] > 0) && ((latticeIndex[c2 - 1] % linearSize) == 0)) {
                latticeIndex[c2 - 1] = 0;
                latticeIndex[c2]++;
                std::cout << "\n";
            }
        }
    }
    std::cout <<"\n****edges****";
    site = 0;
    for (unsigned int ij = 0; ij < numAllNeighbors; ij++) {
        if ((ij % numNearestNeighbors) == 0) { std::cout << "\nsite[" << site++ << "]: "; }
        std::cout << "edge[" << allEdgeNeighbors[ij] << "]->site["<<allSiteNeighbors[ij]<<"] ";
    }
    std::cout << "\n";
}
