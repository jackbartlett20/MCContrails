#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>
#include <atomic>
#include <omp.h>
#include "simulation.h"
#include "population.h"
#include "environment.h"
#include "general.h"
#include "constants.h"

// Runs the simulation
void Simulation::run() {
    read_simulation();
    set_rng(rng_seed_read);
    current_time = 0;
    if (num_writes > 0) {
        write_interval = int_time/num_writes;
    }
    else {
        // next_write_time will be after simulation finishes
        write_interval = int_time;
    }
    next_write_time = write_interval;
    coag_interval = num_dt_for_coag*dt;
    next_coag_time = coag_interval;
    prepare_output();
    
    pop.assign(); // Called after rng set
    
    std::cout << "Starting simulation." << std::endl;
    
    // Start time integration
    while (current_time < int_time) {
        current_time += dt;
        //std::cout << "Current time: " << current_time << std::endl;
        env.set_env(current_time);
        update_water_vol();
        growth();
        if (do_coagulation == 1) {
            if (current_time >= next_coag_time) {
                coagulation();
                next_coag_time += coag_interval;
            }
        }
        freezing();
        if (current_time >= next_write_time) {
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
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> do_coagulation;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> num_dt_for_coag;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> rng_seed_read;
    file.close();
}

// Updates volume of water in each droplet/crystal to reflect change in density
void Simulation::update_water_vol() {
    double dens_ratio;
    // Droplets
    dens_ratio = env.water_density_old / env.water_density;
    #pragma omp parallel for
    for (Superparticle& sp : pop.droplet_sps) {
        double water_vol = dens_ratio * sp.vol * (1 - sp.f_dry);
        double dry_vol = sp.vol * sp.f_dry;
        double new_vol = dry_vol + water_vol;
        sp.vol = new_vol;
        sp.f_dry = dry_vol / new_vol;
    }

    // Crystals
    dens_ratio = env.ice_density_old / env.ice_density;
    #pragma omp parallel for
    for (Superparticle& sp : pop.crystal_sps) {
        double water_vol = dens_ratio * sp.vol * (1 - sp.f_dry);
        double dry_vol = sp.vol * sp.f_dry;
        double new_vol = dry_vol + water_vol;
        sp.vol = new_vol;
        sp.f_dry = dry_vol / new_vol;
    }
}

// Grows the superparticles
void Simulation::growth() {
    // Droplets
    double sum_gn_droplet = 0;
    // Parallelise for loop
    #pragma omp parallel for reduction(+:sum_gn_droplet)
    for (Superparticle& sp : pop.droplet_sps) {
        double growth_rate = growth_rate_liquid(sp.vol, sp.f_dry, sp.kappa);
        sp.f_dry = (sp.f_dry * sp.vol) / (sp.vol + growth_rate * dt);
        sp.vol += growth_rate * dt;
        sum_gn_droplet += growth_rate * sp.n;
        // Check volume
        check_valid_vol(sp.vol);
        // Check f_dry and adjust within tolerance
        sp.f_dry = check_valid_f_dry(sp.f_dry);
    }

    // Crystals
    double sum_gn_crystal = 0;
    // Parallelise for loop
    #pragma omp parallel for reduction(+:sum_gn_crystal)
    for (Superparticle& sp : pop.crystal_sps) {
        double growth_rate = growth_rate_crystal(sp.vol);
        sp.f_dry = (sp.f_dry * sp.vol) / (sp.vol + growth_rate * dt);
        sp.vol += growth_rate * dt;
        sum_gn_crystal += growth_rate * sp.n;
        // Check volume
        check_valid_vol(sp.vol);
        // Check f_dry and adjust within tolerance
        sp.f_dry = check_valid_f_dry(sp.f_dry);
    }

    // Update vapour pressure
    env.Pvap -= (BOLTZMANN_CONSTANT * env.T * (sum_gn_droplet/env.H2O_vol_liquid + sum_gn_crystal/env.H2O_vol_ice)) * dt;
}

// Calculates growth rate for droplets (m3 s-1)
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

    double kelvin_term = std::exp((2*env.sigma_water*env.water_molar_vol)/(IDEAL_GAS_CONSTANT*env.T*r));
    if (std::isinf(kelvin_term)) {
        kelvin_term = std::numeric_limits<double>::max();
    }

    double S_droplet = raoult_term * kelvin_term;

    double diffusivity_mod = env.diffusivity / (r/(r + 0.7*env.mfp_air)
                                                + env.diffusivity/(r*accom_coeff) 
                                                * std::sqrt(2*PI*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T)));

    double k_air_mod = env.k_air / (r/(r + 0.7*env.mfp_air)
                                    + env.k_air/(r*accom_coeff*env.air_density*CP_AIR)
                                    * std::sqrt(2*PI*AIR_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T)));

    double F_d = (env.water_density * IDEAL_GAS_CONSTANT * env.T) / (env.Psat_l * diffusivity_mod * WATER_MOLAR_MASS);

    double F_k = (env.l_v * env.water_density)/(k_air_mod * env.T) * (env.l_v*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.T) - 1);

    // Determine growth rate; set to 0 if in highly volatile regime
    double growth_rate;
    if (raoult_term == 0 && kelvin_term > 1e2) {
        growth_rate = 0;
    }
    else {
        growth_rate = (4*PI*r) * 1/(F_d + F_k) * (env.S_l - S_droplet);
    }
    return growth_rate;
}

// Calculates growth rate for crystals (m3 s-1)
double Simulation::growth_rate_crystal(const double v) {
    const double r = std::pow(3*v/(4*PI), 1./3.);
    const double accom_coeff = 1;
    // Only Kelvin term for crystals
    double S_crystal = std::exp((2*env.sigma_ice*env.ice_molar_vol)/(IDEAL_GAS_CONSTANT*env.T*r));
    double correction_factor = 1 + accom_coeff * env.vapour_thermal_speed * r / (4 * env.diffusivity);
    double J = (PI * std::pow(r,2) * accom_coeff * env.vapour_thermal_speed * env.n_sat) / correction_factor * (env.S_i - S_crystal);
    double growth_rate = env.H2O_vol_ice * J;
    return growth_rate;
}

// Raises error is volume is negative
void Simulation::check_valid_vol(double vol) {
    if (vol <= 0) {
        std::cerr << "Found droplet volume of " << vol << " after growth." << std::endl;
        std::cerr << "Try reducing dt. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Adjusts new f_dry within tolerance; raises error if invalid
double Simulation::check_valid_f_dry(double f_dry) {
    if (f_dry > 1 && f_dry < 1+pop.f_dry_tol) {
        f_dry = 1;
    }
    else if (f_dry < 0 && f_dry > -pop.f_dry_tol) {
        f_dry = 0;
    }
    // Check valid
    if (f_dry < 0 || f_dry > 1) {
        std::cerr << "Error: Found dry fraction of " << f_dry << " after growth." << std::endl;
        std::cerr << "Try reducing dt. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    return f_dry;
}

// Coagulates the superparticles
void Simulation::coagulation() {
    int num_droplet_sps = pop.droplet_sps.size();
    int num_crystal_sps = pop.crystal_sps.size();

    // Initialise coagulation flag vectors (std::atomic<bool> ensures thread-safe updates)
    std::vector<std::atomic<bool>> droplet_sps_coag_flag(num_droplet_sps);
    for (int i = 0; i < num_droplet_sps; i++) {
        droplet_sps_coag_flag[i].store(false);
    }
    std::vector<std::atomic<bool>> crystal_sps_coag_flag(num_crystal_sps);
    for (int i = 0; i < num_crystal_sps; i++) {
        crystal_sps_coag_flag[i].store(false);
    }
    // Initialise empty vectors for created superparticles properties
    std::vector<SPTemp> new_droplet_props;
    std::vector<SPTemp> new_crystal_props;

    // Distribution for choosing random numbers in {0, 1}
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Constant part of probability expression
    // Assumes all superparticles have same number density
    // Which they should for 1-to-1 coagulation
    const double const_fac = pop.n_tot/pop.num_sps * num_dt_for_coag*dt;

    // Droplets with droplets
    // Parallelise outer for loop
    #pragma omp parallel for schedule(dynamic) default(none) \
            shared(num_droplet_sps, distribution, const_fac, pop, droplet_sps_coag_flag, new_droplet_props)
    for (int i = 0; i < num_droplet_sps-1; i++) {
        // Ignore if i has already coagulated
        if (droplet_sps_coag_flag[i].load()) {
            continue;
        }
        for (int j = i+1; j < num_droplet_sps; j++) {
            // Ignore if i or j have already coagulated
            if (droplet_sps_coag_flag[i].load() || droplet_sps_coag_flag[j].load()) {
                continue;
            }
            const Superparticle& sp_i = pop.droplet_sps.at(i);
            const Superparticle& sp_j = pop.droplet_sps.at(j);
            double beta_ij = coag_coeff(sp_i.vol, sp_j.vol);
            double random_number = distribution(global_rng());
            // P_ij = beta_ij * const_fac;
            if (random_number < beta_ij * const_fac) {
                // Coagulation occurs
                droplet_sps_coag_flag[i].store(true);
                droplet_sps_coag_flag[j].store(true);
                SPTemp props = new_props(sp_i, sp_j);
                // Initialise a new droplet superparticle
                #pragma omp critical
                {
                    new_droplet_props.push_back(props);
                }
                break;
            }
        }
    }

    // Crystals with crystals
    // Parallelise outer for loop
    #pragma omp parallel for schedule(dynamic) default(none) \
            shared(num_crystal_sps, distribution, const_fac, pop, crystal_sps_coag_flag, new_crystal_props)
    for (int i = 0; i < num_crystal_sps-1; i++) {
        // Ignore if i has already coagulated
        if (crystal_sps_coag_flag[i].load()) {
            continue;
        }
        for (int j = i+1; j < num_crystal_sps; j++) {
            // Ignore if i or j have already coagulated
            if (crystal_sps_coag_flag[i].load() || crystal_sps_coag_flag[j].load()) {
                continue;
            }
            const Superparticle& sp_i = pop.crystal_sps.at(i);
            const Superparticle& sp_j = pop.crystal_sps.at(j);
            double beta_ij = coag_coeff(sp_i.vol, sp_j.vol);
            double random_number = distribution(global_rng());
            // P_ij = beta_ij * const_fac;
            if (random_number < beta_ij * const_fac) {
                // Coagulation occurs
                crystal_sps_coag_flag[i].store(true);
                crystal_sps_coag_flag[j].store(true);
                SPTemp props = new_props(sp_i, sp_j);
                // Initialise a new crystal superparticle
                #pragma omp critical
                {
                    new_crystal_props.push_back(props);
                }
                break;
            }
        }
    }

    // Droplets with crystals
    // Parallelise outer for loop
    #pragma omp parallel for schedule(dynamic) default(none) \
            shared(num_droplet_sps, num_crystal_sps, distribution, const_fac, pop,\
                   droplet_sps_coag_flag, crystal_sps_coag_flag, new_droplet_props, new_crystal_props)
    for (int i = 0; i < num_droplet_sps; i++) {
        // Ignore if i has already coagulated
        if (droplet_sps_coag_flag[i].load()) {
            continue;
        }
        for (int j = 0; j < num_crystal_sps; j++) {
            // Ignore if i or j have already coagulated
            if (droplet_sps_coag_flag[i].load() || crystal_sps_coag_flag[j].load()) {
                continue;
            }
            const Superparticle& sp_i = pop.droplet_sps.at(i);
            const Superparticle& sp_j = pop.crystal_sps.at(j);
            double beta_ij = coag_coeff(sp_i.vol, sp_j.vol);
            double random_number = distribution(global_rng());
            // P_ij = beta_ij * const_fac;
            if (random_number < beta_ij * const_fac) {
                // Coagulation occurs
                droplet_sps_coag_flag[i].store(true);
                crystal_sps_coag_flag[j].store(true);
                SPTemp props = new_props(sp_i, sp_j);
                // Assume resulting particles are frozen
                // Initialise a new crystal superparticle
                #pragma omp critical
                {
                    new_crystal_props.push_back(props);
                }
                break;
            }
        }
    }

    // Remove coagulated droplet superparticles
    for (int i = num_droplet_sps - 1; i >=0; i--) {
        if (droplet_sps_coag_flag[i].load()) {
            pop.droplet_sps.erase(pop.droplet_sps.begin() + i);
        }
    }
    // Remove coagulated crystal superparticles
    for (int i = num_crystal_sps - 1; i >=0; i--) {
        if (crystal_sps_coag_flag[i].load()) {
            pop.crystal_sps.erase(pop.crystal_sps.begin() + i);
        }
    }
    // Add new droplet superparticles
    for (const SPTemp& props : new_droplet_props) {
        pop.droplet_sps.push_back(Superparticle(pop.sp_ID_count++, props.n, props.vol, props.f_dry, props.kappa, props.ice_germs, false));
    }
    // Add new crystal superparticles
    for (const SPTemp& props : new_crystal_props) {
        pop.crystal_sps.push_back(Superparticle(pop.sp_ID_count++, props.n, props.vol, props.f_dry, props.kappa, props.ice_germs, true));
    }

    // Re-evaluate n_tot and num_sps
    pop.update_n_tot();
    pop.update_num_sps();
}

// Calculates the Brownian coagulation coefficient (m3 s-1) from Fuchs' interpolation formula (Seinfeld and Pandis Table 13.1)
double Simulation::coag_coeff(double vi, double vj) {
    const double ri = std::pow(3*vi/(4*PI), 1./3.);
    const double rj = std::pow(3*vj/(4*PI), 1./3.);
    double Di = diffusivity(ri);
    double Dj = diffusivity(rj);
    double ci = thermal_speed(vi);
    double cj = thermal_speed(vj);
    double gi = coag_g(ri, Di, ci);
    double gj = coag_g(rj, Dj, cj);
    double denom = ((ri + rj)/(ri + rj + std::sqrt(std::pow(gi, 2) + std::pow(gj, 2)))) + ((4*(Di + Dj))/(std::sqrt(std::pow(ci, 2) + std::pow(cj, 2)) * (ri + rj)));
    double beta_ij = 4*PI * (Di+Dj) * (ri + rj) / denom;
    return beta_ij;
}

// Calculates diffusivity (m2 s-1)
double Simulation::diffusivity(double r) {
    double Cc = cscf(r);
    double D = BOLTZMANN_CONSTANT * env.T * Cc / (6 * PI * env.air_viscosity * r);
    return D;
}

// Calculates Cunningham slip correction factor
double Simulation::cscf(double r) {
    double Kn = env.mfp_air / r;
    double Cc = 1 + Kn * (1.257 + 0.4*std::exp(-1.1/Kn));
    return Cc;
}

// Calculates thermal speed of particles with volume v (m s-1); assumes same density as water
double Simulation::thermal_speed(double v) {
    double c = std::sqrt(8 * BOLTZMANN_CONSTANT * env.T / (PI * env.water_density * v));
    return c;
}

// Calculates g for Fuchs coagulation (see Seinfeld and Pandis Table 13.1)
double Simulation::coag_g(double r, double D, double c) {
    double l = 8*D/(PI*c);
    double d = 2*r;
    double g = std::sqrt(2)/(3*d*l) * (std::pow(d + l, 3) - std::pow(std::pow(d, 2) + std::pow(l, 2), 1.5)) - d;
    return g;
}

// Calculates the resulting particle properties from a coagulation and returns them in an SPTemp object
SPTemp Simulation::new_props(Superparticle sp_i, Superparticle sp_j) {
    double new_n = std::min(sp_i.n, sp_j.n);
    double new_vol = sp_i.vol + sp_j.vol;
    double new_f_dry = (sp_i.vol*sp_i.f_dry + sp_j.vol*sp_j.f_dry)/new_vol;
    double new_kappa = (sp_i.vol*sp_i.kappa + sp_j.vol*sp_j.kappa)/new_vol;
    double new_ice_germs = sp_i.ice_germs + sp_j.ice_germs;
    return SPTemp(new_n, new_vol, new_f_dry, new_kappa, new_ice_germs);
}

// Updates number of ice germs in droplets and freezes if >1
void Simulation::freezing() {
    // Number of ice germs formed per droplet vol per second
    double ice_germ_rate = 1e6 * std::exp(-3.5714 * env.T + 858.719);
    // Update number of ice germs
    #pragma omp parallel for
    for (Superparticle& sp : pop.droplet_sps) {
        sp.ice_germs += ice_germ_rate * sp.vol * (1 - sp.f_dry) * dt;
    }
    // Remove if frozen (iterate backwards to allow erase to work)
    double dens_ratio = env.water_density/env.ice_density;
    for (int i = pop.droplet_sps.size() - 1; i >= 0; i--) {
        Superparticle& sp = pop.droplet_sps.at(i);
        if (sp.ice_germs >= 1) {
            // Update properties
            sp.isFrozen = true;
            double water_vol = dens_ratio * sp.vol * (1 - sp.f_dry);
            double dry_vol = sp.vol * sp.f_dry;
            double new_vol = dry_vol + water_vol;
            sp.vol = new_vol;
            sp.f_dry = dry_vol / new_vol;
            // Add to crystal vector
            pop.crystal_sps.push_back(sp);
            // Remove from droplet vector
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
            if (r >= r_output.at(i) && r < r_output.at(i+1)) {
                n_droplet_output.at(i) += sp.n;
                break;
            }
        }
    }

    // Add distribution to n_crystal_output
    for (const Superparticle& sp : pop.crystal_sps) {
        double r = std::pow(3*sp.vol/(4*PI), 1./3.);
        for (int i = 0; i < num_r_intervals_output; i++) {
            if (r >= r_output.at(i) && r < r_output.at(i+1)) {
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