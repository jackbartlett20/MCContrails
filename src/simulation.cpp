#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>
#include <atomic>
#include <omp.h>
#include <yaml-cpp/yaml.h>
#include "simulation.h"
#include "population.h"
#include "environment.h"
#include "general.h"
#include "constants.h"

// Runs the simulation
void Simulation::run(std::string input_path) {
    read_simulation(input_path);
    set_rng(rng_seed_read);
    
    prepare_output();
    env.initialise(input_path);
    pop.assign(input_path, max_sps, num_r_choices); // Uses rng so must be after set_rng

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
    
    std::cout << "Starting simulation." << std::endl;
    
    // Start time integration
    while (current_time < int_time) {
        current_time += dt;
        //std::cout << "Current time: " << current_time << std::endl;
        env.set_env(current_time);
        update_water_vol();
        if (env.get_S_l() >= min_S_l) {
            growth();
        }
        if (do_coagulation == 1) {
            if (current_time >= next_coag_time) {
                coagulation();
                next_coag_time += coag_interval;
            }
        }
        if (env.get_T() < 273.15) {
            freezing();
        }
        if (current_time >= next_write_time) {
            output();
            next_write_time += write_interval;
        }
    }
}

// Reads input file with variables relevant to Simulation
void Simulation::read_simulation(std::string input_path) {
    YAML::Node input_file = YAML::LoadFile(input_path);
    
    int_time               = input_file["simulation"]["int_time"].as<double>();
    dt                     = input_file["simulation"]["dt"].as<double>();
    num_writes             = input_file["simulation"]["num_writes"].as<int>();
    r_output_min           = input_file["simulation"]["r_output_min"].as<double>();
    r_output_max           = input_file["simulation"]["r_output_max"].as<double>();
    num_r_intervals_output = input_file["simulation"]["num_r_intervals_output"].as<int>();
    max_sps                = input_file["simulation"]["max_sps"].as<int>();
    num_r_choices          = input_file["simulation"]["num_r_choices"].as<int>();
    min_S_l                = input_file["simulation"]["min_S_l"].as<double>();
    do_coagulation         = input_file["simulation"]["do_coagulation"].as<int>();
    num_dt_for_coag        = input_file["simulation"]["num_dt_for_coag"].as<int>();
    rng_seed_read          = input_file["simulation"]["rng_seed"].as<unsigned long long>();

    // Check valid
    if (int_time <= 0) {
        std::cerr << "Error: Read in integration time of " << int_time << "." << std::endl;
        std::cerr << "Integration time must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (dt <= 0) {
        std::cerr << "Error: Read in time step of " << dt << "." << std::endl;
        std::cerr << "Time step must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_writes < 0) {
        std::cerr << "Error: Read in " << num_writes << " writes." << std::endl;
        std::cerr << "Number of writes must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_output_min <= 0) {
        std::cerr << "Error: Read in smallest radius in output spectrum of " << r_output_min << "." << std::endl;
        std::cerr << "Smallest radius in output spectrum must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_output_max <= r_output_min) {
        std::cerr << "Error: Read in largest radius in output spectrum of " << r_output_max << "." << std::endl;
        std::cerr << "Largest radius in output spectrum must be greater than smallest radius in output spectrum. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_r_intervals_output <= 0) {
        std::cerr << "Error: Read in " << num_r_intervals_output << " intervals in output spectrum." << std::endl;
        std::cerr << "Number of intervals in output spectrum must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (max_sps <= 0) {
        std::cerr << "Error: Read in total number of superparticles of " << max_sps << "." << std::endl;
        std::cerr << "Total number of superparticles must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_r_choices <= 0) {
        std::cerr << "Error: Read in " << num_r_choices << " radii to choose from during initialisation." << std::endl;
        std::cerr << "Number of radii to choose from must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (min_S_l < 0) {
        std::cerr << "Error: Read in minimum saturation ratio of " << min_S_l << "." << std::endl;
        std::cerr << "Minimum saturation ratio must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (do_coagulation != 0 && do_coagulation != 1) {
        std::cerr << "Error: Read in \"do coagulation?\" of" << do_coagulation << "." << std::endl;
        std::cerr << "Enter 1 for coagulation, 0 for no coagulation. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_dt_for_coag < 0) {
        std::cerr << "Error: Read in number of time steps for coagulation of " << num_dt_for_coag << "." << std::endl;
        std::cerr << "Number of time steps for coagulation must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Updates volume of water in each droplet/crystal to reflect change in density
void Simulation::update_water_vol() {
    double dens_ratio;
    // Droplets
    dens_ratio = env.get_water_density_old() / env.get_water_density();
    #pragma omp parallel for
    for (Superparticle& sp : pop.droplet_sps) {
        sp.vol = sp.dry_vol + dens_ratio * (sp.vol - sp.dry_vol);
    }

    // Crystals
    dens_ratio = env.get_ice_density_old() / env.get_ice_density();
    #pragma omp parallel for
    for (Superparticle& sp : pop.crystal_sps) {
        sp.vol = sp.dry_vol + dens_ratio * (sp.vol - sp.dry_vol);
    }
}

// Grows the superparticles
void Simulation::growth() {
    // Droplets
    double sum_gn_droplet = 0;
    // Parallelise for loop
    #pragma omp parallel for reduction(+:sum_gn_droplet)
    for (Superparticle& sp : pop.droplet_sps) {
        double growth_rate = growth_rate_liquid(sp.vol, sp.dry_vol, sp.kappa);
        sp.vol += growth_rate * dt;
        sum_gn_droplet += growth_rate * sp.n;
        // Check volume
        sp.vol = check_valid_vol(sp.vol, sp.dry_vol);
    }

    // Crystals
    double sum_gn_crystal = 0;
    // Parallelise for loop
    #pragma omp parallel for reduction(+:sum_gn_crystal)
    for (Superparticle& sp : pop.crystal_sps) {
        double growth_rate = growth_rate_crystal(sp.vol);
        sp.vol += growth_rate * dt;
        sum_gn_crystal += growth_rate * sp.n;
        // Check volume
        sp.vol = check_valid_vol(sp.vol, sp.dry_vol);
    }

    // Update vapour pressure
    double delta_Pvap = -(BOLTZMANN_CONSTANT * env.get_T() * (sum_gn_droplet/env.get_H2O_vol_liquid()
                                                              + sum_gn_crystal/env.get_H2O_vol_ice())) * dt;
    env.update_Pvap_after_growth(delta_Pvap);
}

// Calculates growth rate for droplets (m3 s-1)
double Simulation::growth_rate_liquid(const double v, const double v_dry, const double kappa) {
    const double r = v_to_r(v);
    const double f_dry = v_dry/v;
    const double accom_coeff = 1;

    double raoult_term;
    if (kappa == 0) {
        raoult_term = 1;
    }
    else {
        raoult_term = (1 - f_dry)/(1 - (1 - kappa)*f_dry);
    }

    double kelvin_term = std::exp((2*env.get_sigma_water()*env.get_water_molar_vol())/(IDEAL_GAS_CONSTANT*env.get_T()*r));
    if (std::isinf(kelvin_term)) {
        kelvin_term = std::numeric_limits<double>::max();
    }

    double S_droplet = raoult_term * kelvin_term;

    double diffusivity_mod = env.get_diffusivity() / (r/(r + 0.7*env.get_mfp_air())
                                                + env.get_diffusivity()/(r*accom_coeff) 
                                                * std::sqrt(2*PI*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.get_T())));

    double k_air_mod = env.get_k_air() / (r/(r + 0.7*env.get_mfp_air())
                                    + env.get_k_air()/(r*accom_coeff*env.get_air_density()*CP_AIR)
                                    * std::sqrt(2*PI*AIR_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.get_T())));

    double F_d = (env.get_water_density() * IDEAL_GAS_CONSTANT * env.get_T()) / (env.get_Psat_l() * diffusivity_mod * WATER_MOLAR_MASS);

    double F_k = (env.get_l_v() * env.get_water_density())/(k_air_mod * env.get_T()) * (env.get_l_v()*WATER_MOLAR_MASS/(IDEAL_GAS_CONSTANT*env.get_T()) - 1);

    double growth_rate = (4*PI*r) * 1/(F_d + F_k) * (env.get_S_l() - S_droplet);
    return growth_rate;
}

// Calculates growth rate for crystals (m3 s-1)
double Simulation::growth_rate_crystal(const double v) {
    const double r = v_to_r(v);
    const double accom_coeff = 1;
    // Only Kelvin term for crystals
    double S_crystal = std::exp((2*env.get_sigma_ice()*env.get_ice_molar_vol())/(IDEAL_GAS_CONSTANT*env.get_T()*r));
    double correction_factor = 1 + accom_coeff * env.get_vapour_thermal_speed() * r / (4 * env.get_diffusivity());
    double J = (PI * std::pow(r,2) * accom_coeff * env.get_vapour_thermal_speed() * env.get_n_sat()) / correction_factor * (env.get_S_i() - S_crystal);
    double growth_rate = env.get_H2O_vol_ice() * J;
    return growth_rate;
}

// Raises error is volume is negative
double Simulation::check_valid_vol(double vol, double dry_vol) {
    if (vol < dry_vol) {
        //std::cerr << "Found volume = " << vol/dry_vol << " * dry volume after growth." << std::endl;
        vol = dry_vol;
    }
    return vol;
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
        pop.droplet_sps.push_back(Superparticle(pop.sp_ID_count++, props.n, props.vol, props.dry_vol, props.kappa, props.ice_germs, false));
    }
    // Add new crystal superparticles
    for (const SPTemp& props : new_crystal_props) {
        pop.crystal_sps.push_back(Superparticle(pop.sp_ID_count++, props.n, props.vol, props.dry_vol, props.kappa, props.ice_germs, true));
    }

    // Re-evaluate n_tot and num_sps
    pop.update_n_tot();
    pop.update_num_sps();
}

// Calculates the Brownian coagulation coefficient (m3 s-1) from Fuchs' interpolation formula (Seinfeld and Pandis Table 13.1)
double Simulation::coag_coeff(double vi, double vj) {
    const double ri = v_to_r(vi);
    const double rj = v_to_r(vj);
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
    double D = BOLTZMANN_CONSTANT * env.get_T() * Cc / (6 * PI * env.get_air_viscosity() * r);
    return D;
}

// Calculates Cunningham slip correction factor
double Simulation::cscf(double r) {
    double Kn = env.get_mfp_air() / r;
    double Cc = 1 + Kn * (1.257 + 0.4*std::exp(-1.1/Kn));
    return Cc;
}

// Calculates thermal speed of particles with volume v (m s-1); assumes same density as water
double Simulation::thermal_speed(double v) {
    double c = std::sqrt(8 * BOLTZMANN_CONSTANT * env.get_T() / (PI * env.get_water_density() * v));
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
    double new_dry_vol = sp_i.dry_vol + sp_j.dry_vol;
    double new_kappa = (sp_i.vol*sp_i.kappa + sp_j.vol*sp_j.kappa)/new_vol;
    double new_ice_germs = sp_i.ice_germs + sp_j.ice_germs;
    return SPTemp(new_n, new_vol, new_dry_vol, new_kappa, new_ice_germs);
}

// Updates number of ice germs in droplets and freezes if >1
void Simulation::freezing() {
    // Update number of ice germs using Koop et al. 2000
    double u_w_i_minus_u_w_0 = 210368. + 131.438 * env.get_T() - 3.32373e6/env.get_T() - 41729.1 * std::log(env.get_T());
    double a_w_i = std::exp(u_w_i_minus_u_w_0/(IDEAL_GAS_CONSTANT * env.get_T()));
    double v_w_0 = -230.76 - 0.1478 * env.get_T() + 4099.2/env.get_T() + 48.8341 * std::log(env.get_T());
    double v_i_0 = 19.43 - 2.2e-3 * env.get_T() + 1.08e-5 * std::pow(env.get_T(), 2);
    #pragma omp parallel for
    for (Superparticle& sp : pop.droplet_sps) {
        //sp.ice_germs += ice_germ_rate * (sp.vol - sp.dry_vol) * dt;
        double r = v_to_r(sp.vol);
        double P_droplet = (env.get_P_ambient() + (2*env.get_sigma_water()/r)) / 1e9; // GPa
        double a_w = (sp.vol - sp.dry_vol)/(sp.vol - (1-sp.kappa)*sp.dry_vol);
        double v_w_minus_v_i = v_w_0 * (P_droplet - 0.5 * 1.6 * std::pow(P_droplet, 2) - 1./6. * -8.8 * std::pow(P_droplet, 3)) 
                               - v_i_0 * (P_droplet - 0.5 * 0.22 * std::pow(P_droplet, 2) - 1./6. * -0.17 * std::pow(P_droplet, 3));
        double delta_a_w = a_w * std::exp(v_w_minus_v_i/(IDEAL_GAS_CONSTANT * env.get_T())) - a_w_i;
        double J;
        if (delta_a_w < 0.26 || delta_a_w > 0.34) {
            J = 0;
        }
        else {
            J = 1e6 * std::pow(10, (-906.7 + 8502 * delta_a_w - 26924. * std::pow(delta_a_w, 2) + 29180. * std::pow(delta_a_w, 3)));
        }
        sp.ice_germs += J * (sp.vol - sp.dry_vol) * dt;
    }
    // Remove if frozen (iterate backwards to allow erase to work)
    double dens_ratio = env.get_water_density()/env.get_ice_density();
    for (int i = pop.droplet_sps.size() - 1; i >= 0; i--) {
        Superparticle& sp = pop.droplet_sps.at(i);
        if (sp.ice_germs >= 1) {
            // Update properties
            sp.isFrozen = true;
            sp.vol = sp.dry_vol + dens_ratio * (sp.vol - sp.dry_vol);
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

    // Create output directory if needed
    std::filesystem::path outputDir = "./output";
    if (!std::filesystem::exists(outputDir)) {
        if (std::filesystem::create_directories(outputDir)) {
                std::cout << "Created output directory." << std::endl;
        }
        else {
            std::cerr << "Failed to create output directory at '" << outputDir.string() << "'. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    std::cout << "Output will be written to '" << outputDir.string() << "'." << std::endl;
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
        double r = v_to_r(sp.vol);
        for (int i = 0; i < num_r_intervals_output; i++) {
            if (r >= r_output.at(i) && r < r_output.at(i+1)) {
                n_droplet_output.at(i) += sp.n;
                break;
            }
        }
    }

    // Add distribution to n_crystal_output
    for (const Superparticle& sp : pop.crystal_sps) {
        double r = v_to_r(sp.vol);
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

    std::ofstream file; // Variable used for both outputs

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
        << n_crystal_output.at(i) << ", " << dndlogr_crystal_output.at(i) << "\n";
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
    file << current_time << ", " << env.get_T() << ", " << env.get_Pvap() << ", " << env.get_Psat_l() << ", " << env.get_Psat_i() << "\n";
    file.close();

    first_write = false;
}