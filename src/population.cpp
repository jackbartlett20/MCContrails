#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <random>
#include "population.h"
#include "general.h"
#include "constants.h"

void Population::assign() {
    read_population();

    // Read species
    std::vector<Species> species_vec = read_species();
    int num_species = species_vec.size();
    std::cout << "Number of species read: " << num_species << std::endl;

    // Calculate total particle density
    n_tot = 0;
    for (const Species& species : species_vec) {
        n_tot += species.n;
    }
    double n_per_sp = n_tot/max_sps;

    // Assign superparticles
    num_sps = 0;
    sp_ID_count = 0;
    int species_ID = 0;
    for (const Species& species : species_vec) {
        int num_sps_for_species = std::round(species.n/n_per_sp);
        species_ID++;
        std::cout << "Number of superparticles created from species " << species_ID << ": " << num_sps_for_species << std::endl;
        if (num_sps_for_species == 0) {
            continue;
        }
        double logr_mean = std::log(species.r_mean);
        double logr_min = logr_mean - 4*species.SD;
        double logr_max = logr_mean + 4*species.SD;
        std::vector<double> logr_vec = linspace(logr_min, logr_max, num_r_choices);
        std::vector<double> weights = normal_dist(logr_vec, logr_mean, species.SD);
        std::vector<double> vs = choose_vs(num_sps_for_species, logr_vec, weights);
        for (int i = 0; i < num_sps_for_species; i++) {
            droplet_sps.push_back(Superparticle(sp_ID_count++, n_per_sp, vs.at(i), species.f_dry, species.kappa, 0, false));
            num_sps++;
        }
    }
}

// Reads input file with variables relevant to Population
void Population::read_population() {
    std::ifstream file("input/population.in");
    if (!file.is_open()) {
        std::cerr << "Error opening population.in" << std::endl;
        exit(EXIT_FAILURE);
    }
    file >> max_sps;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> num_r_choices;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> f_dry_tol;
    file.close();
}

// Reads each species from input file
std::vector<Species> Population::read_species() {
    std::vector<Species> species_vec;

    std::ifstream file("input/species.in");
    if (!file.is_open()) {
        std::cerr << "Error opening psr.in" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line=="END") {
            break;
        }
        // Read each line segment into string
        std::stringstream ss(line);
        std::string n_str, r_mean_str, SD_str, f_dry_str, kappa_str;
        std::getline(ss, n_str, ',');
        std::getline(ss, r_mean_str, ',');
        std::getline(ss, SD_str, ',');
        std::getline(ss, f_dry_str, ',');
        std::getline(ss, kappa_str, ',');

        // Convert strings to doubles
        double n = std::stod(n_str);
        double r_mean = std::stod(r_mean_str);
        double SD = std::stod(SD_str);
        double f_dry = std::stod(f_dry_str);
        double kappa = std::stod(kappa_str);

        // Check valid
        if (n < 0) {
            std::cerr << "Error: Read in number density of " << n << "." << std::endl;
            std::cerr << "Number density should be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (r_mean <= 0) {
            std::cerr << "Error: Read in mean radius of " << r_mean << "." << std::endl;
            std::cerr << "Mean radius must be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (f_dry < 0 ||  f_dry > 1) {
            std::cerr << "Error: Read in dry fraction of " << f_dry << "." << std::endl;
            std::cerr << "Dry fraction must be between 0 and 1. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (kappa < 0) {
            std::cerr << "Error: Read in hygroscopicity of " << kappa << "." << std::endl;
            std::cerr << "Hygroscopicity must be >= 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }

        species_vec.push_back(Species(n, r_mean, SD, f_dry, kappa));
    }
    file.close();
    return species_vec;
}

// Randomly chooses a set of volumes for superparticles according to weights
std::vector<double> Population::choose_vs(int num_sps_for_species, std::vector<double> logr_range, std::vector<double> weights) {
    std::vector<double> vs(num_sps_for_species);
    //std::random_device rd;
    //std::mt19937 gen(rd()); // Seed Mersenne Twister engine
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    for (int i = 0; i < num_sps_for_species; i++) {
        //int random_index = distribution(gen);
        int random_index = distribution(global_rng());
        double r = std::exp(logr_range.at(random_index));
        vs.at(i) = 4*PI/3*std::pow(r, 3);
    }
    return vs;
}

Species::Species(double n, double r_mean, double SD, double f_dry, double kappa) {
    this->n = n;
    this->r_mean = r_mean;
    this->SD = SD;
    this->f_dry = f_dry;
    this->kappa = kappa;
}

Superparticle::Superparticle(int ID, double n, double vol, double f_dry, double kappa, double ice_germs, bool isFrozen) {
    this->ID = ID;
    this->n = n;
    this->vol = vol;
    this->f_dry = f_dry;
    this->kappa = kappa;
    this->ice_germs = ice_germs;
    this->isFrozen = isFrozen;
}