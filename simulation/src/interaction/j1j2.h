#ifndef J1J2_H
#define J1J2_H

#include <vector>
#include <iomanip>
#include <fstream>
#include "../random_numbers/random_numbers.h"
#include "../lattice/lattice.h"

class J1J2{
private:
	int disorderSeed{-1}; //disorder seed
	int disorderSample{-1};
	unsigned int nrBonds; //nr. bonds
	double p2; //occupation probability of j2 bonds
	double j2; //strenght of the second bond type. The strength of the other bond is J_1=1
	std::vector<unsigned char> bonds; //bond array
    
public:
	J1J2(const Lattice& inLattice, double inP2, const double inJ2, const int inDisorderSeed, const int inDisorderSample);
	~J1J2(){};	

	const unsigned char * get_bonds_data()const{return bonds.data();}
	unsigned char * get_bonds_data(){return bonds.data();}
	const double get_j2()const{return j2;}

	void save_bond_map_for_2d_partition_function(const Lattice& inLattice)const;
	void save_bond_map_for_2d_alternating_xy(const Lattice& inLattice)const;

	void set_j1j2_disorder(int inDisorderSeed, int inDisorderSample, double inP2);

	void set_j1j2_disorder_sample(const unsigned int inDisorderSample);

	void print_2d(const Lattice& inLattice)const;
	void print(const Lattice& inLattice)const;

};

#endif
