#ifndef LATTICE_H
#define LATTICE_H

#include<vector>
#include<cmath>
#include<iostream>

class Lattice{

private:

unsigned short dimension;//dimension of hypercubic lattice
unsigned int linearSize;//linear size of hypercubic lattice
unsigned short numNearestNeighbors;//number of nearest neighbors of a site
unsigned int numSites;//number of all sites of the lattice
unsigned int numEdges;//number of all edges of the lattice
unsigned int numAllNeighbors;//number of all neighboring sites of all sites
unsigned int numAllAdjacentSites;//number of all adjacent sites of all edges
unsigned int numAllCoordinates;//number of all coordinates of all sites

std::vector<unsigned int> allSiteNeighbors;//indexes of all neighboring sites of all sites
std::vector<unsigned int> allEdgeNeighbors;//indexes of all neighboring edges of all edges 
std::vector<unsigned int> allAdjacentSites;//indexes of all adjacent sites of all edges
std::vector<unsigned int> allCoordinates;//coordinates of all sites 
std::vector<unsigned int> dimensionFactors;//factors of dimension for calculating neighboring sites and coordinates



public:
Lattice();
Lattice(const unsigned short inDimension, const unsigned int inLinearSize);
~Lattice();

void print()const;

unsigned short get_dimension()const{return dimension;}
unsigned int get_linear_size()const{return linearSize;}
unsigned short get_num_nearest_neighbors()const {return numNearestNeighbors;}
unsigned int get_num_sites()const{return numSites;}
unsigned int get_num_edges()const{return numEdges;}
unsigned int get_num_all_neighbors()const{return numAllNeighbors;}
unsigned int get_num_all_adjacent_sites()const{return numAllAdjacentSites;}
unsigned int get_num_all_coordinates()const{return numAllCoordinates;}

const unsigned int * get_all_site_neighbors_data()const{return allSiteNeighbors.data();}
const unsigned int * get_all_edge_neighbors_data()const{return allEdgeNeighbors.data();}
const unsigned int * get_all_adjacent_sites_data()const{return allAdjacentSites.data();}
const unsigned int * get_all_coordinates_data()const{return allCoordinates.data();}
const unsigned int * get_dimension_factors_data()const{return dimensionFactors.data();}


};

#endif
