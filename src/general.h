#ifndef GENERAL
#define GENERAL

#include <vector>

std::vector<double> linspace(double start, double stop, int n);

std::vector<double> expspace(double start, double stop, int n);

std::vector<double> normal_dist(std::vector<double> x, double x_mean, double sigma);

#endif