#ifndef ISINGCONFIGURATION_H
#define ISINGCONFIGURATION_H


#include <vector>
#include "../lattice/lattice.h"
#include "../random_numbers/random_numbers.h"



class IsingConfiguration{

private:

unsigned int numSpins;
std::vector<signed char> spins;


public:

IsingConfiguration();

~IsingConfiguration(){}

IsingConfiguration(const Lattice& inLattice); 
void print(const Lattice & inLattice)const;
void set_spins(const int inS=1);
void set_spins_rand(const double inP, RandomNumbers& inRandomNumbers);
unsigned int get_num_spins()const{ return numSpins; }

const signed char * get_spins_data()const{return spins.data();}
signed char * get_spins_data(){return spins.data();}


};
#endif
