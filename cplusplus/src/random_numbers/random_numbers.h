#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

#include <iostream>
#include <vector>
#include "../lattice/lattice.h"

std::vector<unsigned char> bernoulli_vector(double inP, unsigned int inSize, int inDisorderSeed, int inDisorderSample);



#ifdef USE_MKL


#include <mkl_vsl.h>

class RandomNumbers {
private:
    unsigned int size;
    unsigned int index;
    double* randomBuffer01;
    double* randomDoubles;
    int* randomIntegers;
    VSLStreamStatePtr generator;

public:
    RandomNumbers(unsigned int inSize, int inSeed, int inThreadNum = 0) {
        size = inSize;
        index = size;
        randomBuffer01 = new double[size];
        randomDoubles  = new double[size];
        randomIntegers = new int[size];
        for (unsigned int c = 0; c < size; c++) {
            randomBuffer01[c] = -1.0;
            randomDoubles[c]  = -1.0;
            randomIntegers[c] = -1;
        }
        vslNewStream(&generator, VSL_BRNG_MT2203 + inThreadNum, inSeed);
    }

    RandomNumbers(const Lattice& inLattice, int inSeed, int inThreadNum = 0) {
        size = 4 * ((inLattice.get_num_edges() < 100000) ? 100000 : inLattice.get_num_edges());
        index = size;
        randomBuffer01 = new double[size];
        randomDoubles  = new double[size];
        randomIntegers = new int[size];
        for (unsigned int c = 0; c < size; c++) {
            randomBuffer01[c] = -1.0;
            randomDoubles[c]  = -1.0;
            randomIntegers[c] = -1;
        }
        vslNewStream(&generator, VSL_BRNG_MT2203 + inThreadNum, inSeed);
    }

    ~RandomNumbers() {
        delete[] randomBuffer01;
        delete[] randomDoubles;
        delete[] randomIntegers;
        vslDeleteStream(&generator);
    }

    unsigned int get_index() { return index; }
    double* get_random_buffer01_data() { return randomBuffer01; }
    double* get_random_doubles_data()  { return randomDoubles; }
    int*    get_random_integers_data() { return randomIntegers; }

    void set_index(const unsigned int inCurrentIndex) {
        if (inCurrentIndex > size) { std::cout << "\nIndex too large"; }
        index = inCurrentIndex;
    }

    void fill_random_buffer_01(unsigned int& currentIndex, const unsigned int numRequestedNumbers) {
        if ((size - currentIndex) < numRequestedNumbers) {
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, generator, currentIndex, randomBuffer01, 0.0, 1.0);
            currentIndex = 0;
        }
    }

    void generate_uniform_random_integers(const unsigned int inNumRequestedNumbers, int inMinInteger, int inMaxInteger) {
        if (inNumRequestedNumbers <= size) {
            viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, generator, inNumRequestedNumbers, randomIntegers, inMinInteger, inMaxInteger + 1);
        } else {
            std::cout << "\nRequested number of random integers exceeds the size of memory.";
        }
    }

    void print_random_buffer() {
        std::cout << "\nRandom Buffer:\n";
        for (unsigned int c = 0; c < size; c++)
            std::cout << "randomBuffer01[" << c << "]=" << randomBuffer01[c] << "\n";
    }

    std::vector<unsigned char> bernoulli_vector(double inP, unsigned int inSize);
};


#else  // USE_MKL not defined — use standard <random>

#include <random>

class RandomNumbers {
private:
    unsigned int size;
    unsigned int index;
    double* randomBuffer01;
    double* randomDoubles;
    int* randomIntegers;
    std::mt19937 generator;

public:
    RandomNumbers(unsigned int inSize, int inSeed, int inThreadNum = 0) {
        size = inSize;
        index = size;
        randomBuffer01 = new double[size];
        randomDoubles  = new double[size];
        randomIntegers = new int[size];
        for (unsigned int c = 0; c < size; c++) {
            randomBuffer01[c] = -1.0;
            randomDoubles[c]  = -1.0;
            randomIntegers[c] = -1;
        }
        generator.seed(inSeed + inThreadNum);
    }

    RandomNumbers(const Lattice& inLattice, int inSeed, int inThreadNum = 0) {
        size = 4 * ((inLattice.get_num_edges() < 100000) ? 100000 : inLattice.get_num_edges());
        index = size;
        randomBuffer01 = new double[size];
        randomDoubles  = new double[size];
        randomIntegers = new int[size];
        for (unsigned int c = 0; c < size; c++) {
            randomBuffer01[c] = -1.0;
            randomDoubles[c]  = -1.0;
            randomIntegers[c] = -1;
        }
        generator.seed(inSeed + inThreadNum);
    }

    ~RandomNumbers() {
        delete[] randomBuffer01;
        delete[] randomDoubles;
        delete[] randomIntegers;
    }

    unsigned int get_index() { return index; }
    double* get_random_buffer01_data() { return randomBuffer01; }
    double* get_random_doubles_data()  { return randomDoubles; }
    int*    get_random_integers_data() { return randomIntegers; }

    void set_index(const unsigned int inCurrentIndex) {
        if (inCurrentIndex > size) { std::cout << "\nIndex too large"; }
        index = inCurrentIndex;
    }

    void fill_random_buffer_01(unsigned int& currentIndex, const unsigned int numRequestedNumbers) {
        if ((size - currentIndex) < numRequestedNumbers) {
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            for (unsigned int c = 0; c < currentIndex; c++)
                randomBuffer01[c] = dist(generator);
            currentIndex = 0;
        }
    }

    void generate_uniform_random_integers(const unsigned int inNumRequestedNumbers, int inMinInteger, int inMaxInteger) {
        if (inNumRequestedNumbers <= size) {
            std::uniform_int_distribution<int> dist(inMinInteger, inMaxInteger);
            for (unsigned int c = 0; c < inNumRequestedNumbers; c++)
                randomIntegers[c] = dist(generator);
        } else {
            std::cout << "\nRequested number of random integers exceeds the size of memory.";
        }
    }

    void print_random_buffer() {
        std::cout << "\nRandom Buffer:\n";
        for (unsigned int c = 0; c < size; c++)
            std::cout << "randomBuffer01[" << c << "]=" << randomBuffer01[c] << "\n";
    }

    std::vector<unsigned char> bernoulli_vector(double inP, unsigned int inSize);
};

#endif // USE_MKL

#endif // RANDOMNUMBERS_H
