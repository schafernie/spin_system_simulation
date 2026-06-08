#ifndef CLUSTERS_H
#define CLUSTERS_H

#include "../lattice/lattice.h"

class Clusters{

private:

unsigned int numSites;
unsigned int numEdges;
std::vector<int> clusterIds;
std::vector<int> stack;
std::vector<unsigned char> occupiedEdges;


public:

~Clusters(){}
Clusters(const Lattice & inLattice)
{
    numSites = inLattice.get_num_sites();
    numEdges = inLattice.get_num_edges();
    clusterIds.resize(numSites);
    stack.resize(numSites);
    occupiedEdges.resize(numEdges);
}


void set_cluster_ids(const int inValue=-1){
    for(int & e : clusterIds ){
         e = inValue;
    }
}

const unsigned int get_num_sites()const{return numSites;}
const unsigned int get_num_edges()const{return numEdges;}

int * get_stack_data(){return stack.data();}
int * get_cluster_ids_data(){return clusterIds.data();}
unsigned char * get_occupied_edges_data(){return occupiedEdges.data();}


};
#endif

