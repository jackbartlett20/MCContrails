#ifndef GENERAL
#define GENERAL

#include <vector>
#include <random>

static unsigned long long rng_seed;

std::vector<double> linspace(double start, double stop, int n);

std::vector<double> expspace(double start, double stop, int n);

std::vector<double> normal_dist(std::vector<double> x, double x_mean, double sigma);

std::mt19937_64& global_rng();

#endif