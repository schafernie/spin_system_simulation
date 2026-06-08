#ifndef STATISTICS_H
#define STATISTICS_H

#include <cmath>
#include <vector>

// Arithmetic mean
template<typename T>
double mean(const std::vector<T>& inValues) {
    double result{ 0.0 };
    for (const T& e : inValues)
        result += static_cast<double>(e);
    return result / static_cast<double>(inValues.size());
}

// Sample variance.  inDdof=1 gives the unbiased estimator, inDdof=0 the biased one.
template<typename T>
double variance(const std::vector<T>& inValues, const unsigned int inDdof = 1) {
    const double n = static_cast<double>(inValues.size());
    double m1{ 0.0 }, m2{ 0.0 };
    for (const T& e : inValues) {
        const double x = static_cast<double>(e);
        m1 += x;
        m2 += x * x;
    }
    m1 /= n; m2 /= n;
    return (m2 - m1 * m1) * (n / (n - inDdof));
}

// N-th raw moment
template<typename T>
double moment(const std::vector<T>& inValues, const unsigned int nthMoment = 1) {
    double result{ 0.0 };
    for (const T& e : inValues)
        result += std::pow(static_cast<double>(e), nthMoment);
    return result / static_cast<double>(inValues.size());
}

#endif 
