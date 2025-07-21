#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>
#include "simulation.h"
#include "population.h"
#include "environment.h"
#include "general.h"
#include "constants.h"

// Runs the simulation
void Simulation::run() {
    read_simulation();
    current_time = 0;
    if (num_writes > 0) {
        write_interval = int_time/num_writes;
    }
    else {
        // next_write_time will be after simulation finishes
        write_interval = int_time;
    }
    next_write_time = write_interval;
    prepare_output();
    
    std::cout << "Starting simulation." << std::endl;
    while (current_time < int_time) {
        current_time += dt;
        //std::cout << "Current time: " << current_time << std::endl;
        env.set_env(current_time);
        growth();
        freeze();
        if (current_time>next_write_time) {
            output();
            next_write_time += write_interval;
        }
    }
}

// Reads input file with variables relevant to Simulation
void Simulation::read_simulation() {
    std::ifstream file("input/simulation.in");
    if (!file.is_open()) {
        std::cerr << "Error opening simulation.in" << std::endl;
        exit(EXIT_FAILURE);
    }
    file >> int_time;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> dt;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> num_writes;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> r_output_min;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> r_output_max;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> num_r_intervals_output;
    file.close();
}

// Grows the superparticles
void Simulation::growth() {
    double sum_gn = 0;

    // Droplets
    for (Superparticle& sp : pop.droplet_sps) {
        double growth_rate = growth_rate_liquid(sp.vol, sp.f_dry, sp.kappa);
        sp.f_dry = (sp.f_dry * sp.vol) / (sp.vol + growth_rate * dt);
        sp.vol += growth_rate * dt;
        sum_gn += growth_rate * sp.n;
        
        // Adjust f_dry within tolerance
        if (sp.f_dry > 1 && sp.f_dry < 1+pop.f_dry_tol) {
            sp.f_dry = 1;
        }
        else if (sp.f_dry < 0 && sp.f_dry > -pop.f_dry_tol) {
            sp.f_dry = 0;
        }
        // Check valid
        if (sp.f_dry < 0 || sp.f_dry > 1) {
            std::cerr << "Error: Found dry fraction of " << sp.f_dry << " after growth." << std::endl;
            std::cerr << "Try reducing dt. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Crystals
    for (Superparticle& sp : pop.crystal_sps) {
        double growth_rate = growth_rate_crystal(sp.vol);
        sp.f_dry = (sp.f_dry * sp.vol) / (sp.vol + growth_rate * dt);
        sp.vol += growth_rate * dt;
        sum_gn += growth_rate * sp.n;
        
        // Adjust f_dry within tolerance
        if (sp.f_dry > 1 && sp.f_dry < 1+pop.f_dry_tol) {
            sp.f_dry = 1;
        }
        else if (sp.f_dry < 0 && sp.f_dry > -pop.f_dry_tol) {
            sp.f_dry = 0;
        }
        // Check valid
        if (sp.f_dry < 0 || sp.f_dry > 1) {
            std::cerr << "Error: Found dry fraction of " << sp.f_dry << " after growth." << std::endl;
            std::cerr << "Try reducing dt. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Update vapour pressure
    env.Pvap += (BOLTZMANN_CONSTANT * env.T * (-sum_gn / WATER_MOLECULAR_VOLUME)) * dt;
}

// Calculates growth rate for droplets
double Simulation::growth_rate_liquid(const double v, const double f_dry, const double kappa) {
    const double r = std::pow(3*v/(4*PI), 1./3.);
    const double accom_coeff = 1;
    double raoult_term;
    if (kappa == 0) {
        raoult_term = 1;
    }
    else {
        raoult_term = (1 - f_dry)/(1 - (1 - kappa)*f_dry);
    } 

    double kelvin_term = std::exp((2*env.sigma_water*WATER_MOLAR_MASS)/(IDEAL_GAS_CONSTANT*env.T*WATER_DENSITY*r));

    double S_droplet = raoult_term * kelvin_term;

    double diffusivity_mod = env.diffusivity / (r/(r + 0.7*env.mfp_air)
                                                + env.diffusivity/(r*accom_coeff) 
                                                * std::sqrt(2*PI*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T)));

    double k_air_mod = env.k_air / (r/(r + 0.7*env.mfp_air)
                                    + env.k_air/(r*accom_coeff*env.air_density*CP_AIR)
                                    * std::sqrt(2*PI*AIR_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T)));

    double F_d = (WATER_DENSITY * IDEAL_GAS_CONSTANT * env.T) / (env.Psat_l * diffusivity_mod * WATER_MOLAR_MASS);

    double F_k = (env.l_v * WATER_DENSITY)/(k_air_mod * env.T) * (env.l_v*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T) - 1);

    double growth_rate = (4*PI*r) * 1/(F_d + F_k) * (env.S_l - S_droplet);
    return growth_rate;
}

// Calculates growth rate for crystals
double Simulation::growth_rate_crystal(const double v) {
    const double r = std::pow(3*v/(4*PI), 1./3.);
    const double accom_coeff = 1;
    double correction_factor = 1 + accom_coeff * env.vapour_thermal_speed * r / (4 * env.diffusivity);
    double J = (PI * std::pow(r,2) * accom_coeff * env.vapour_thermal_speed * (env.S_i - 1) * env.n_sat) / correction_factor;
    double growth_rate = WATER_MOLECULAR_VOLUME * J;
    return growth_rate;
}

// Updates number of ice germs in droplets and freezes if >1
void Simulation::freeze() {
    // Number of ice germs formed per droplet vol per second
    double ice_germ_rate = 1e6 * std::exp(-3.5714 * env.T + 858.719);
    for (int i = pop.droplet_sps.size() - 1; i >= 0; i--) {
        Superparticle& sp = pop.droplet_sps.at(i);
        sp.ice_germs += ice_germ_rate * sp.vol * (1 - sp.f_dry) * dt;
        if (sp.ice_germs >= 1) {
            sp.isFrozen = true;
            pop.crystal_sps.push_back(sp);
            pop.droplet_sps.erase(pop.droplet_sps.begin() + i);
        }
    }
}

// Prepares variables which are used at each output time
void Simulation::prepare_output() {
    first_write = true;
    r_output.resize(num_r_intervals_output+1);
    r_m_output.resize(num_r_intervals_output);
    n_droplet_output.resize(num_r_intervals_output);
    dndlogr_droplet_output.resize(num_r_intervals_output);
    n_crystal_output.resize(num_r_intervals_output);
    dndlogr_crystal_output.resize(num_r_intervals_output);

    // Find boundaries
    r_output = expspace(r_output_min, r_output_max, num_r_intervals_output+1);

    // Find midpoints
    for (int i = 0; i < num_r_intervals_output; i++) {
        r_m_output.at(i) = 0.5 * (r_output.at(i) + r_output.at(i+1));
    }

    // Calculate dlogr
    // Since r_output is linear in logspace, all dlogr are the same
    dlogr_output = std::log(r_output.at(1)) - std::log(r_output.at(0));
}

// Outputs PSDs and environmental variables to relevant files
void Simulation::output() {
    // Set n and dndlogr to zero
    for (int i = 0; i < num_r_intervals_output; i++) {
        n_droplet_output.at(i) = 0;
        dndlogr_droplet_output.at(i) = 0;
        n_crystal_output.at(i) = 0;
        dndlogr_crystal_output.at(i) = 0;
    }

    // Add distribution to n_droplet_output
    for (const Superparticle& sp : pop.droplet_sps) {
        double r = std::pow(3*sp.vol/(4*PI), 1./3.);
        for (int i = 0; i < num_r_intervals_output; i++) {
            if (r > r_output.at(i) && r < r_output.at(i+1)) {
                n_droplet_output.at(i) += sp.n;
                break;
            }
        }
    }

    // Add distribution to n_crystal_output
    for (const Superparticle& sp : pop.crystal_sps) {
        double r = std::pow(3*sp.vol/(4*PI), 1./3.);
        for (int i = 0; i < num_r_intervals_output; i++) {
            if (r > r_output.at(i) && r < r_output.at(i+1)) {
                n_crystal_output.at(i) += sp.n;
                break;
            }
        }
    }

    // Calculate dn/dlogr
    for (int i = 0; i < num_r_intervals_output; i++) {
        dndlogr_droplet_output.at(i) = n_droplet_output.at(i) / dlogr_output;
        dndlogr_crystal_output.at(i) = n_crystal_output.at(i) / dlogr_output;
    }

    std::ofstream file; // Variable used for all outputs

    // Write PSDs
    if (first_write == true) {
        file.open("output/psd.out");
    }
    else {
        file.open("output/psd.out", std::ios::app);
    }
    if (!file.is_open()) {
        std::cerr << "Error opening psd.out" << std::endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < num_r_intervals_output; i++) {
        file << current_time << ", " << r_m_output.at(i) << ", " << n_droplet_output.at(i) << ", " << dndlogr_droplet_output.at(i) << ", "
        << n_crystal_output.at(i) << ", " << dndlogr_crystal_output.at(i) << std::endl;
    }
    file.close();

    // Write environment variables
    if (first_write == true) {
        file.open("output/environment.out");
    }
    else {
        file.open("output/environment.out", std::ios::app);
    }
    if (!file.is_open()) {
        std::cerr << "Error opening environment.out" << std::endl;
        exit(EXIT_FAILURE);
    }
    file << current_time << ", " << env.T << ", " << env.Pvap << ", " << env.Psat_l << ", " << env.Psat_i << std::endl;
    file.close();

    first_write = false;
}