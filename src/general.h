#ifndef GENERAL
#define GENERAL

#include <vector>
#include <random>

std::vector<double> linspace(double start, double stop, int n);

std::vector<double> expspace(double start, double stop, int n);

std::vector<double> normal_dist(std::vector<double> x, double x_mean, double sigma);

double r_to_v(double r);
double v_to_r(double v);

void set_rng(unsigned long long rng_seed_read);

std::mt19937_64& global_rng();

#endif