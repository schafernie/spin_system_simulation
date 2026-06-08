#include <iostream>
#include <string>

#include "../simulation/simulation.h"

using namespace std;

int main(int argc, char* argv[]) {
    
    Parameters inParameters;
    int option{0}; 
    bool printHelp{false};
    unsigned int c{1};

    while ((c < argc) && (argv[c][0] == '-')) {
        string arg = argv[c++];
        if (arg == "-L")
            inParameters.linearSize = static_cast<unsigned int>(stoul(argv[c]));
        else if (arg == "-D")
            inParameters.dimension = static_cast<short>(stoi(argv[c]));
        else if (arg == "-O")
            option = static_cast<short>(stoi(argv[c]));
        else if (arg == "-T")
            inParameters.temperature = static_cast<double>(stod(argv[c]));
        else if (arg == "-P")
            inParameters.p2 = static_cast<double>(stod(argv[c]));
        else if (arg == "-J")
            inParameters.j2 = static_cast<double>(stod(argv[c]));
        else if (arg == "--DisorderSeed")
            inParameters.disorderSeed = static_cast<int>(stoi(argv[c]));
        else if (arg == "--DisorderSample")
            inParameters.disorderSample = static_cast<int>(stoi(argv[c]));
        else if (arg == "--ThermalSeed")
            inParameters.thermalSeed = static_cast<int>(stoi(argv[c]));
        else if (arg == "--EquilibrationTime")
            inParameters.thermalEquilibrationTime = static_cast<unsigned int>(stoul(argv[c]));
        else if (arg == "--ThermalSamples")
            inParameters.thermalSamples = static_cast<unsigned int>(stoul(argv[c]));
        else if (arg == "--TimeBetweenSamples")
            inParameters.thermalTimeBetweenSamples = static_cast<unsigned int>(stoul(argv[c]));
        else if (arg == "--NrDisorderSamples")
            inParameters.nrDisorderSamples = static_cast<unsigned int>(stoul(argv[c]));
        else if (arg == "--FixedDisorder") {
            inParameters.fixedDisorder = true;
            c--;  // no value follows this flag
        }
        else {
            cerr << " unknown option:" << arg << "\n";
            printHelp = true;
            break;
        }
        c++;
    }

    if (printHelp || (argc - c) != 0) {
        cerr << "USAGE: " << argv[0] << " {<options>}\n"
             << "  -O <int>                 option (0=equilibration, 1=disorder sampling, -1=debugging)\n"
             << "  -D <int>                 dimension\n"
             << "  -L <int>                 linear size\n"
             << "  -T <double>              temperature\n"
             << "  -P <double>              bond occupation probability\n"
             << "  -J <double>              j2 coupling\n"
             << "  --DisorderSeed <int>     disorder seed\n"
             << "  --DisorderSample <int>   disorder sample\n"
             << "  --ThermalSeed <int>         thermal seed\n"
             << "  --EquilibrationTime <int>   equilibration time\n"
             << "  --ThermalSamples <int>      thermal samples per disorder realization\n"
             << "  --TimeBetweenSamples <int>  MC steps between samples\n"
             << "  --NrDisorderSamples <int>   number of disorder realizations\n"
             << "  --FixedDisorder             keep disorder sample fixed across all iterations\n";
        return 1;
    }
    if(inParameters.thermalSeed == 0)
        inParameters.thermalSeed = inParameters.disorderSeed+inParameters.disorderSample;  
    if (option == 0){
        j1j2_equilibration(inParameters);
    }
    else if (option == 1) {
 
        j1j2_disorder_sampling(inParameters);
    }
    else if (option == -1) {
        debugging(inParameters);
    }

    return 0;
}

