#include "j1j2.h"


J1J2::J1J2(const Lattice& inLattice, double inP2, const double inJ2, const int inDisorderSeed, const int inDisorderSample)
{                       
    nrBonds = inLattice.get_num_edges();
    disorderSeed=inDisorderSeed;
    disorderSample=inDisorderSample;
    p2 = inP2;
    j2 = inJ2;
    bonds.resize(nrBonds);
    set_j1j2_disorder(disorderSeed, disorderSample, p2);
}


void J1J2::set_j1j2_disorder(int inDisorderSeed, int inDisorderSample, double inP2){
    std::vector<unsigned char> bernoullis = bernoulli_vector(inP2, nrBonds, inDisorderSeed, inDisorderSample);
    for (unsigned int c = 0; c < nrBonds; c++) {
        if(bernoullis[c]){bonds[c] = 0;}
        else { bonds[c] = 1;}
    }
}

void J1J2::set_j1j2_disorder_sample(const unsigned int inDisorderSample){
     std::vector<unsigned char> bernoullis = bernoulli_vector(p2, nrBonds, disorderSeed, inDisorderSample);
    for (unsigned int c = 0; c < nrBonds; c++) {
        if(bernoullis[c]){bonds[c] = 0;}
        else { bonds[c] = 1;}
    }
}




void J1J2::print(const Lattice& inLattice)const {
    const unsigned short dimension{ inLattice.get_dimension() };
    const unsigned linearSize{inLattice.get_linear_size() };
    const unsigned short numNearestNeighbors = inLattice.get_num_nearest_neighbors();
    const unsigned int  numAllNeighbors = inLattice.get_num_all_neighbors();
    const unsigned int * allSiteNeighbors = inLattice.get_all_site_neighbors_data();
    const unsigned int* allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();

    std::cout << "\nj1--j2 disordered ferromagnet: "; 
    std::cout << "\nj2: " << j2;
    std::cout << "\np2: " << p2;
    std::cout << "\ndimension: " << dimension;
    std::cout << "\nlinear size: " << linearSize;
    std::cout << "\ndisorder seed: " << disorderSeed;
    std::cout << "\ndisorder sample: " << disorderSample;
    std::cout << "\n1 = j1 bond, 0 = j2 bond \n";


    unsigned int site{ 0 };
    for (unsigned int ij = 0; ij < numAllNeighbors; ij++) {
        if ((ij % numNearestNeighbors) == 0) { std::cout << "\nsite[" << site++ << "]: "; }
        std::cout << "bond[" << allEdgeNeighbors[ij] << "] = " << static_cast<int>(bonds[allEdgeNeighbors[ij]]) << " "
        <<" to site "<<allSiteNeighbors[ij]<<";  ";
    }
    std::cout << "\n";
}


void J1J2::print_2d(const Lattice& inLattice)const {
    const unsigned short dimension{ inLattice.get_dimension() };
    const unsigned linearSize{inLattice.get_linear_size() };
    const unsigned short numNearestNeighbors = inLattice.get_num_nearest_neighbors();
    const unsigned int  numAllNeighbors = inLattice.get_num_all_neighbors();
    const unsigned int * allSiteNeighbors = inLattice.get_all_site_neighbors_data();
    const unsigned int* allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();

    if(dimension != 2){
        std::cout << "\nprint_2d only implemented for 2d lattices!";
        return;
    }
    std::cout << "\nj1--j2 disordered ferromagnet: "; 
    std::cout << "\np2: " << j2;
    std::cout << "\np2: " << p2;
    std::cout << "\ndimension: " << dimension;
    std::cout << "\nlinear size: " << linearSize;
    std::cout << "\ndisorder seed: " << disorderSeed;
    std::cout << "\ndisorder sample: " << disorderSample;
    std::cout << "\n1 = j1 bond, 0 = j2 bond \n";

    unsigned int site{ 0 };
    for (unsigned int ij = 0; ij < numAllNeighbors; ij++) {
        if ((ij % numNearestNeighbors) == 0) { std::cout << "\nsite[" << site++ << "]: "; }
        std::cout << "bond[" << allEdgeNeighbors[ij] << "] = " << static_cast<int>(bonds[allEdgeNeighbors[ij]]) << " "
        <<" to site "<<allSiteNeighbors[ij]<<";  ";
    }

    
    
    int counter, offset;
    std::cout << "\n\nindices of bonds where even bonds are horizontal and odd bonds are vertical:\n";
    for (int y=0; y < linearSize; y++){
        offset=y*linearSize;
        counter=offset;
        for (int x0=0; x0 < linearSize; x0++){
        std::cout << std::setw(6) <<offset+counter;
        counter+=2; 
        }
        std::cout<<std::endl;
        counter=2*y*linearSize;
        for (int x1=0; x1 < linearSize; x1++){
        std::cout << std::setw(6) <<counter+1; 
        counter+=2;
        }
         std::cout<<std::endl;
    }
    std::cout << "\ncorresponding bonds:\n";
    for (int y=0; y < linearSize; y++){
        offset=y*linearSize;
        counter=offset;
        for (int x0=0; x0 < linearSize; x0++){
        std::cout << std::setw(6) << static_cast<int>(bonds[offset+counter]);
        counter+=2; 
        }
        std::cout<<std::endl;
        counter=2*y*linearSize;
        for (int x1=0; x1 < linearSize; x1++){
        std::cout << std::setw(6) << static_cast<int>(bonds[counter+1]); 
        counter+=2;
        }
         std::cout<<std::endl;
    }
}



void J1J2::save_bond_map_for_2d_partition_function(const Lattice& inLattice)const {
    const unsigned short dimension{ inLattice.get_dimension() };
    const unsigned linearSize{ inLattice.get_linear_size() };
    const unsigned int* allEdgeNeighbors = inLattice.get_all_edge_neighbors_data();

    if(dimension != 2){
        std::cout << "\nsave_bond_map_for_2d_partition_function only implemented for 2d lattices!";
        return;
    }


    std::stringstream stringStream;
    stringStream <<"J1J2BondZ" << "_D" << dimension << "_L" << linearSize << "_P" <<
    std::fixed << std::setprecision(4) << p2 <<"_J"<< j2
    <<"_DisorderSeed"<<disorderSeed<<"_DisorderSample"<<disorderSample<<"_.dat";
    std::ofstream file("./../data/bonds/"+ stringStream.str());
    file << std::setprecision(15) <<  std::defaultfloat;
    file << linearSize <<"    "<<linearSize<< "\n";
    unsigned int counter{0};
    for (unsigned int y = 0; y < linearSize; y++) {
        for (unsigned int x = 0; x < linearSize; x++) {
            file <<std::setw(8)<< x <<std::setw(8)<< y <<std::setw(8)<<"E";
            if(bonds[allEdgeNeighbors[counter]]){file<<std::setw(20)<<1;}
            else{file<<std::setw(20) << std::fixed <<j2;}
            file<<"\n";
            counter+=2;
            file <<std::setw(8)<< x <<std::setw(8)<< y <<std::setw(8)<<"S"<<std::setw(8);
            if(bonds[allEdgeNeighbors[counter]]){file<<std::setw(20)<<1;}
            else{file<<std::setw(20) << std::fixed <<j2;}
            file<<"\n";
            counter+=2;

        }
    }
    file.close(); 
}

void J1J2::save_bond_map_for_2d_alternating_xy(const Lattice& inLattice)const {
    const unsigned short dimension{ inLattice.get_dimension() };
    const unsigned linearSize{ inLattice.get_linear_size() };

    if(dimension != 2){
        std::cout << "\nsave_bond_map_for_2d_partition_function only implemented for 2d lattices!";
        return;
    }

     std::stringstream stringStream;

    stringStream <<"J1J2Bond" << "_D" << dimension << "_L" << linearSize << "_P" <<
    std::fixed << std::setprecision(4) << p2 <<"_J"<< j2 <<"_DisorderSeed"<<disorderSeed<<"_DisorderSample"<<disorderSample<<"_.dat";
    std::ofstream file("./../data/bonds/"+ stringStream.str());
    
    file <<"#Alternating rows describing horizontal and vertical bonds, respectively, starting with the horizontal bonds of the first row."
         <<"\n#1 = j1 bond, 0 = j2 bond\n";

    unsigned int counter{0};
    unsigned int offset;
    std::cout << "\n";
    for (int y=0; y < linearSize; y++){
        offset=y*linearSize;
        counter=offset;
        for (int x0=0; x0 < linearSize; x0++){
        file<<static_cast<int>(bonds[offset+counter])<<" ";
        counter+=2; 
        }
        file<<std::endl;
        counter=2*y*linearSize;
        for (int x1=0; x1 < linearSize; x1++){
        file<<static_cast<int>(bonds[counter+1])<<" "; 
        counter+=2;
        }
         file<<std::endl;
    }
    file.close();
}