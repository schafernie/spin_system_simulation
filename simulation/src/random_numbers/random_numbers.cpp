#include "random_numbers.h"


#ifdef USE_MKL


std::vector<unsigned char> bernoulli_vector(double inP, unsigned int inSize, int inDisorderSeed, int inDisorderSample) {
    std::vector<unsigned char> bernoullis(inSize);
    VSLStreamStatePtr randomGenerator;
    vslNewStream(&randomGenerator, VSL_BRNG_SFMT19937, inDisorderSeed);
    vslSkipAheadStream(randomGenerator, (MKL_LONG)inSize * inDisorderSample);
    std::vector<int> randomNumbers(inSize);
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, randomGenerator, inSize, randomNumbers.data(), inP);
    for (unsigned int c = 0; c < inSize; c++)
        bernoullis[c] = static_cast<unsigned char>(randomNumbers[c]);
    vslDeleteStream(&randomGenerator);
    return bernoullis;
}

std::vector<unsigned char> RandomNumbers::bernoulli_vector(double inP, unsigned int inSize) {
    std::vector<unsigned char> bernoullis(inSize);
    std::vector<int> randomNumbers(inSize);
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, generator, inSize, randomNumbers.data(), inP);
    for (unsigned int c = 0; c < inSize; c++)
        bernoullis[c] = static_cast<unsigned char>(randomNumbers[c]);
    return bernoullis;
}


#else  // USE_MKL not defined — use standard <random>


std::vector<unsigned char> bernoulli_vector(double inP, unsigned int inSize, int inDisorderSeed, int inDisorderSample) {
    std::vector<unsigned char> bernoullis(inSize);
    std::mt19937 generator(inDisorderSeed);
    std::bernoulli_distribution dist(inP);
    for (unsigned long c = 0; c < (unsigned long)inSize * inDisorderSample; c++)
        dist(generator);
    for (unsigned int c = 0; c < inSize; c++)
        bernoullis[c] = dist(generator) ? 1 : 0;
    return bernoullis;
}

std::vector<unsigned char> RandomNumbers::bernoulli_vector(double inP, unsigned int inSize) {
    std::vector<unsigned char> bernoullis(inSize);
    std::bernoulli_distribution dist(inP);
    for (unsigned int c = 0; c < inSize; c++)
        bernoullis[c] = dist(generator) ? 1 : 0;
    return bernoullis;
}

#endif // USE_MKL
