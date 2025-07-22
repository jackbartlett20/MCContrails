#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "general.h"
#include "constants.h"

// Returns a vector of n linearly spaced values inclusive of start and stop
std::vector<double> linspace(double start, double stop, int n) {
    double step = (stop - start)/(n-1);
    std::vector<double> vec(n);
    for (int i = 0; i < n; i++) {
        vec.at(i) = start + step*i;
    }
    return vec;
}

// Returns a vector of n exponentially spaced values inclusive of start and stop
std::vector<double> expspace(double start, double stop, int n) {
    double log_start = std::log(start);
    double log_stop = std::log(stop);
    double log_step = (log_stop - log_start)/(n-1);
    std::vector<double> vec(n);
    for (int i = 0; i < n; i++) {
        vec.at(i) = std::exp(log_start + log_step*i);
    }
    return vec;
}

// Calculates a normal distribution from the input vector, returns vector
std::vector<double> normal_dist(std::vector<double> x, double x_mean, double sigma) {
    int vec_size = x.size();
    std::vector<double> normal(vec_size);
    for (int i = 0; i < vec_size; i++) {
        normal.at(i) = 1/(std::sqrt(2*PI)*sigma) * std::exp(-std::pow(x.at(i)-x_mean, 2)/(2*std::pow(sigma, 2)));
    }
    return normal;
}

// Returns a reference to the global random number generator; only initialised when this function is first called
std::mt19937_64& global_rng() {
    // rng_seed is set in Simulation::set_rng()
    static std::mt19937_64 rng(rng_seed);
    return rng;
}