#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <random>
#include <yaml-cpp/yaml.h>
#include "population.h"
#include "general.h"
#include "constants.h"

void Population::assign(std::string input_path, int max_sps, int num_r_choices) {
    this->max_sps = max_sps;
    this->num_r_choices = num_r_choices;

    // Read species
    std::vector<Species> species_vec = read_species(input_path);
    int num_species = species_vec.size();
    std::cout << "Number of species read: " << num_species << std::endl;

    // Calculate total particle density
    n_tot = 0;
    for (const Species& species : species_vec) {
        n_tot += species.n;
    }
    double n_per_sp = n_tot/max_sps;

    // Assign superparticles
    sp_ID_count = 0;
    int species_ID = 0;
    for (const Species& species : species_vec) {
        int num_sps_for_species = std::round(species.n/n_per_sp);
        species_ID++;
        std::cout << "Number of superparticles created from species " << species_ID << ": " << num_sps_for_species << std::endl;
        if (num_sps_for_species == 0) {
            continue;
        }
        double logr_mean = std::log(species.GMR);
        double SD = std::log(species.GSD);
        double logr_min = logr_mean - 4*SD;
        double logr_max = logr_mean + 4*SD;
        std::vector<double> logr_vec = linspace(logr_min, logr_max, num_r_choices);
        std::vector<double> weights = normal_dist(logr_vec, logr_mean, SD);
        std::vector<double> vs = choose_vs(num_sps_for_species, logr_vec, weights);
        for (int i = 0; i < num_sps_for_species; i++) {
            droplet_sps.push_back(Superparticle(sp_ID_count++, n_per_sp, vs.at(i), vs.at(i)*species.f_dry, species.kappa, 0, false));
        }
    }
    // Determine num_sps
    update_num_sps();
}

// Reads each species from input file
std::vector<Species> Population::read_species(std::string input_path) {
    std::vector<Species> species_vec;

    YAML::Node input_file = YAML::LoadFile(input_path);

    // Iterate over species
    for (const auto& speciesNode : input_file["species"]) {
        double n     = speciesNode["n"].as<double>();
        double GMR   = speciesNode["GMR"].as<double>();
        double GSD   = speciesNode["GSD"].as<double>();
        double f_dry = speciesNode["f_dry"].as<double>();
        double kappa = speciesNode["kappa"].as<double>();

        // Check valid
        if (n < 0) {
            std::cerr << "Error: Read in number density of " << n << "." << std::endl;
            std::cerr << "Number density should be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (GMR <= 0) {
            std::cerr << "Error: Read in geometric mean radius of " << GMR << "." << std::endl;
            std::cerr << "Geometric mean radius must be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (f_dry < 0 || f_dry > 1) {
            std::cerr << "Error: Read in dry fraction of " << f_dry << "." << std::endl;
            std::cerr << "Dry fraction must be between 0 and 1. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (kappa < 0) {
            std::cerr << "Error: Read in hygroscopicity of " << kappa << "." << std::endl;
            std::cerr << "Hygroscopicity must be >= 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }

        species_vec.push_back(Species(n, GMR, GSD, f_dry, kappa));
    }
    return species_vec;
}

// Randomly chooses a set of volumes for superparticles according to weights
std::vector<double> Population::choose_vs(int num_sps_for_species, std::vector<double> logr_range, std::vector<double> weights) {
    std::vector<double> vs(num_sps_for_species);
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    for (int i = 0; i < num_sps_for_species; i++) {
        int random_index = distribution(global_rng());
        double r = std::exp(logr_range.at(random_index));
        vs.at(i) = r_to_v(r);
    }
    return vs;
}

// Sets n_tot inside Population
void Population::update_n_tot() {
    double sum = 0;
    for (Superparticle& sp : droplet_sps) {
        sum += sp.n;
    }
    for (Superparticle& sp : crystal_sps) {
        sum += sp.n;
    }
    n_tot = sum;
}

// Sets num_sps inside Population
void Population::update_num_sps() {
    num_sps = droplet_sps.size() + crystal_sps.size();
}

Species::Species(double n, double GMR, double GSD, double f_dry, double kappa) {
    this->n = n;
    this->GMR = GMR;
    this->GSD = GSD;
    this->f_dry = f_dry;
    this->kappa = kappa;
}

Superparticle::Superparticle(int ID, double n, double vol, double dry_vol, double kappa, double ice_germs, bool isFrozen) {
    this->ID = ID;
    this->n = n;
    this->vol = vol;
    this->dry_vol = dry_vol;
    this->kappa = kappa;
    this->ice_germs = ice_germs;
    this->isFrozen = isFrozen;
}

SPTemp::SPTemp(double n, double vol, double dry_vol, double kappa, double ice_germs) {
    this->n = n;
    this->vol = vol;
    this->dry_vol = dry_vol;
    this->kappa = kappa;
    this->ice_germs = ice_germs;
}