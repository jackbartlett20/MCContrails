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
#include "params.h"

void Population::assign(Params& params, const double initial_n_to_EI_ratio) {
    max_sps = params.max_sps;
    num_r_choices = params.num_r_choices;

    std::vector<Species> species_vec = params.species_vec;
    int num_species = species_vec.size();
    std::cout << "Number of species read: " << num_species << std::endl;

    // Calculate total initial particle density
    n_tot = 0;
    for (const Species& species : species_vec) {
        n_tot += species.EI * initial_n_to_EI_ratio;
    }
    double n_per_sp = n_tot/max_sps;

    // Assign superparticles
    sp_ID_count = 0;
    int species_ID = 0;
    for (const Species& species : species_vec) {
        int num_sps_for_species = std::round(species.EI*initial_n_to_EI_ratio/n_per_sp);
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
            sps.push_back(Superparticle(sp_ID_count++, n_per_sp, vs.at(i), vs.at(i)*species.f_dry, species.kappa, 0, false));
        }
    }
    // Determine num_sps
    update_num_sps();
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
    for (Superparticle& sp : sps) {
        sum += sp.n;
    }
    n_tot = sum;
}

// Sets num_sps inside Population
void Population::update_num_sps() {
    num_sps = sps.size();
}

Species::Species(double EI, double GMR, double GSD, double f_dry, double kappa) {
    this->EI = EI;
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