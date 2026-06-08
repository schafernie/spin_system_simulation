#include "simulation.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>





void j1j2_equilibration(const Parameters& inParameters) {
    const unsigned short dimension = inParameters.dimension;
    const unsigned int linearSize = inParameters.linearSize;
    const double temperature = inParameters.temperature;
    const double beta = 1.0 / temperature;
    const double p2 = inParameters.p2;
    const double j2 = inParameters.j2;
    const int disorderSeed = inParameters.disorderSeed;
    const int disorderSample = inParameters.disorderSample;
    const int thermalSeed = inParameters.thermalSeed;
    const unsigned int thermalEquilibrationTime = inParameters.thermalEquilibrationTime;

    cout << "\nj1j2_equilibration()";
    cout << "\ndimension: " << dimension;
    cout << "\nlinearSize: " << linearSize;
    cout << "\ntemperature: " << temperature;
    cout << "\np2: " << p2;
    cout << "\nj2: " << j2;
    cout << "\ndisorderSeed: " << disorderSeed;
    cout << "\ndisorderSample: " << disorderSample;
    cout << "\nthermalEquilibrationTime: " << thermalEquilibrationTime;
    cout <<endl;

    const string path{"./../data/new/"};
    Lattice lattice(dimension, linearSize);
    J1J2 j1j2(lattice, p2, j2, disorderSeed, disorderSample);
    RandomNumbers randomNumbers(lattice, thermalSeed);
    Clusters clusters(lattice);
    IsingConfiguration configuration(lattice);

    stringstream ss;
    ss << setprecision(6) << fixed
       << "J1J2Equilibration_D" << dimension << "_L" << linearSize
       << "_P" << fixed << setprecision(4) << p2
       << "_J" << fixed << setprecision(4) << j2
       << "_T" << fixed << setprecision(6) << temperature
       << "_DisorderSeed" <<disorderSeed<<"_DisorderSample"<<disorderSample<< "_ThermalSeed" << thermalSeed;

    // Random start
    configuration.set_spins_rand(0.5, randomNumbers);
    ofstream file(path + ss.str() + "_RS_.dat");
    file << "# temperature: " << temperature << " , disorder seed: " << disorderSeed<<" , disorder sample: " << disorderSample << " , thermal seed: " << thermalSeed << "\n";
    file << "# time , energy , magnetization\n";
    for (unsigned int t = 0; t < thermalEquilibrationTime; t++) {
        file<< t<<" "<< energy(lattice, j1j2, configuration) << " " << magnetization(configuration) << "\n";
        mcmc_swendsen_wang(lattice, j1j2, beta, configuration, clusters, randomNumbers, 1);
    }
    file.close();

    // Ground state start
    configuration.set_spins(1);
    file.open(path + ss.str() + "_GS_.dat");
    file << "# temperature: " << temperature << " , disorder seed: " << disorderSeed <<" , disorder sample: " << disorderSample <<  " , thermal seed: " << thermalSeed << "\n";
    file << "# energy magnetization\n";
    for (unsigned int t = 0; t < thermalEquilibrationTime; t++) {
        file<< t<<" " << energy(lattice, j1j2, configuration) << " " << magnetization(configuration) << "\n";
        mcmc_swendsen_wang(lattice, j1j2, beta, configuration, clusters, randomNumbers, 1);
    }
    file.close();
    cout<<"ending.\n";
}


void j1j2_disorder_sampling(const Parameters& inParameters) {
    const double timeStart = omp_get_wtime();
    const unsigned short dimension = inParameters.dimension;
    const unsigned int linearSize = inParameters.linearSize;
    const double temperature = inParameters.temperature;
    const double beta = 1.0 / temperature;
    const double p2 = inParameters.p2;
    const double j2 = inParameters.j2;
    const int disorderSeed = inParameters.disorderSeed;
    const int disorderSample = inParameters.disorderSample;
    const int thermalSeed = inParameters.thermalSeed;
    const unsigned int thermalEquilibrationTime = inParameters.thermalEquilibrationTime;
    const unsigned int thermalSamples = inParameters.thermalSamples;
    const unsigned int thermalTimeBetweenSamples = inParameters.thermalTimeBetweenSamples;
    const unsigned int nrDisorderSamples = inParameters.nrDisorderSamples;
    const bool fixedDisorder = inParameters.fixedDisorder;

    cout << "\nj1j2_disorder_sampling()";
    cout << "\ndimension: " << dimension;
    cout << "\nlinearSize: " << linearSize;
    cout << "\ntemperature: " << temperature;
    cout << "\np2: " << p2;
    cout << "\nj2: " << j2;
    cout << "\ndisorderSeed: " << disorderSeed;
    cout << "\ndisorderSample: " << disorderSample;
    cout << "\nthermalEquilibrationTime: " << thermalEquilibrationTime;
    cout << "\nthermalSamples: " << thermalSamples;
    cout << "\nthermalTimeBetweenSamples: " << thermalTimeBetweenSamples;
    cout << "\nnrDisorderSamples: " << nrDisorderSamples;

    const string path{"./../data/new/"};
    Lattice lattice(dimension, linearSize);
    const unsigned int numSites = lattice.get_num_sites();
    J1J2 j1j2(lattice, p2, j2, disorderSeed, disorderSample);
    RandomNumbers randomNumbers(lattice, thermalSeed, 0);
    Clusters clusters(lattice);
    IsingConfiguration configuration(lattice);
    configuration.set_spins_rand(0.5, randomNumbers);

    stringstream ss;
    ss << (fixedDisorder ? "J1J2DisorderSamplingFixed" : "J1J2DisorderSampling")
       << "_D" << dimension << "_L" << linearSize
       << "_P" << fixed << setprecision(4) << p2
       << "_J" << fixed << setprecision(4) << j2
       << "_T" << fixed << setprecision(6) << temperature
       << "_ThermalSamples" << thermalSamples
       << "_DisorderSeed" << disorderSeed << "_DisorderSample" << disorderSample << "_.dat";
    cout << "\n" << path + ss.str();

    ofstream file(path + ss.str());
    file << setprecision(15) << defaultfloat;
    file << "# temperature: " << temperature << " , p2: " << p2
         << " , equilibration time: " << thermalEquilibrationTime
         << " , thermal samples: " << thermalSamples
         << " , time between samples: " << thermalTimeBetweenSamples << "\n";
    file << "# disorder sample , |m| , |m|^2 , |m|^3 , |m|^4 , e , e^2 , e^3 , e^4\n";

    std::vector<double> magnetizationSamples(thermalSamples);
    std::vector<double> energySamples(thermalSamples);

    for (unsigned int c = 0; c < nrDisorderSamples; c++) {
        if (!fixedDisorder){j1j2.set_j1j2_disorder_sample(disorderSample + c);}
        mcmc_swendsen_wang(lattice, j1j2, beta, configuration, clusters, randomNumbers, thermalEquilibrationTime);
        for (unsigned int t = 0; t < thermalSamples; t++) {
            mcmc_swendsen_wang(lattice, j1j2, beta, configuration, clusters, randomNumbers, thermalTimeBetweenSamples);
            magnetizationSamples[t] = std::abs(magnetization(configuration)) / static_cast<double>(numSites);
            energySamples[t] = energy(lattice, j1j2, configuration) / static_cast<double>(numSites);
        }
        file << (fixedDisorder ? disorderSample : disorderSample + c) << " "
             << moment(magnetizationSamples, 1) << " "
             << moment(magnetizationSamples, 2) << " "
             << moment(magnetizationSamples, 3) << " "
             << moment(magnetizationSamples, 4) << " "
             << moment(energySamples, 1) << " "
             << moment(energySamples, 2) << " "
             << moment(energySamples, 3) << " "
             << moment(energySamples, 4) << "\n";
    }
    file.close();
    cout << "\nprogram took " << (omp_get_wtime() - timeStart) / 60.0 << " minutes.\n";
    //j1j2.save_bond_map_for_2d_partition_function(lattice);
}





void debugging(const Parameters& inParameters){
    unsigned short dimension = inParameters.dimension;
    unsigned int linearSize = inParameters.linearSize;
    double temperature = inParameters.temperature;
    double beta=1./temperature;
    double p2 = inParameters.p2;
    double j2 = inParameters.j2;
    int disorderSeed = inParameters.disorderSeed;
    int disorderSample = inParameters.disorderSample;

    
    
    
    // Print the parameters
    cout << "\nj1j2_debugging()";
    cout << "\nPrint the parameters";
    cout << "\ndimension: " << dimension;
    cout << "\nlinearSize: " << linearSize;
    cout << "\nbeta: " << beta;
    cout << "\ntemperature: " << temperature;
    cout << "\np2: " << p2;
    cout << "\nj2: " << j2;
    cout << "\ndisorderSeed: " << disorderSeed; 
    cout << "\ndisorderSample: " << disorderSample;


    //Print the lattice
    cout << "\nPrint the lattice";
    Lattice lattice(dimension, linearSize);
    lattice.print();
     
    //Print random numbers
      cout << "\nPrint random numbers";
    RandomNumbers randomNumbers(10, 0, 0);
    randomNumbers.print_random_buffer();
    unsigned int currentIndex = randomNumbers.get_index();
    cout << "\nCurrentIndex: " << currentIndex;
    cout << "\nFilling 1 number ...";
    randomNumbers.fill_random_buffer_01(currentIndex, 1);
    cout << "\nCurrentIndex after filling: " << currentIndex;
    randomNumbers.print_random_buffer();
    double  * randomBuffer01 = randomNumbers.get_random_buffer01_data();
    cout << "\nFirst 9 random numbers:\n";  
    for (unsigned int c = 0; c < 9; c++) {
        cout << "randomBuffer01[" << currentIndex << "]=" << randomBuffer01[currentIndex++] << "\n";
    }
    cout << "\nCurrentIndex: " << currentIndex;
    randomNumbers.print_random_buffer();
    cout << "\nFilling 5 numbers ...";
    randomNumbers.fill_random_buffer_01(currentIndex, 5);
    cout << "\nCurrentIndex after filling 5 numbers: " << currentIndex;
    randomNumbers.print_random_buffer();
    cout << "\nFilling 7 numbers ... ";
    randomNumbers.fill_random_buffer_01(currentIndex, 7);
    randomNumbers.print_random_buffer();
    cout << "\nCurrentIndex after filling 7 numbers: " << currentIndex;
    randomNumbers.print_random_buffer();
    
    // Print the bonds
    cout << "\nPrint the interactions";
    J1J2 j1j2(lattice, p2, j2, disorderSeed, disorderSample);
    j1j2.print_2d(lattice);
    j1j2.save_bond_map_for_2d_partition_function(lattice);
    j1j2.save_bond_map_for_2d_alternating_xy(lattice);

    RandomNumbers randomNumbers2(lattice, 1, 0);
    IsingConfiguration configuration(lattice);
    configuration.set_spins(1);
    configuration.print(lattice);
  
    cout<<"magnetization: "<< magnetization(configuration);
    cout<<"\nenergy: "<< energy(lattice, j1j2, configuration);
    
    configuration.set_spins(-1);
    configuration.print(lattice);
    cout<<"magnetization: "<< magnetization(configuration);
    cout<<"\nenergy: "<< energy(lattice, j1j2, configuration);
    configuration.set_spins_rand(0.5, randomNumbers2);
    configuration.print(lattice);
    cout<<"magnetization: "<< magnetization(configuration);
    cout<<"\nenergy: "<< energy(lattice, j1j2, configuration)<<"\n";

    mcmc_single_spin_flip(lattice, j1j2, beta, configuration, randomNumbers2, 1000);
    configuration.print(lattice);
    cout<<"magnetization: "<< magnetization(configuration);
    cout<<"\nenergy: "<< energy(lattice, j1j2, configuration)<<"\n";

    Clusters clusters(lattice);
    mcmc_swendsen_wang(lattice, j1j2, beta, configuration, clusters, randomNumbers2,1000); 

    configuration.print(lattice);
    cout<<"magnetization: "<< magnetization(configuration);
    cout<<"\nenergy: "<< energy(lattice, j1j2, configuration)<<"\n";


}


