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
    double water_dens_ratio = env.get_water_density_old() / env.get_water_density();
    double ice_dens_ratio = env.get_ice_density_old() / env.get_ice_density();
    #pragma omp parallel for
    for (Superparticle& sp : pop.sps) {
        double dens_ratio = sp.isFrozen ? ice_dens_ratio : water_dens_ratio;
        sp.vol = sp.dry_vol + dens_ratio * (sp.vol - sp.dry_vol);
    }
}

// Grows the superparticles
void Simulation::growth() {
    double sum_gn_droplet = 0;
    double sum_gn_crystal = 0;
    // Parallelise for loop
    #pragma omp parallel for reduction(+:sum_gn_droplet, sum_gn_crystal)
    for (Superparticle& sp : pop.sps) {
        double growth_rate;
        if (sp.isFrozen) {
            growth_rate = growth_rate_crystal(sp.vol);
            sum_gn_crystal += growth_rate * sp.n;
        }
        else {
            growth_rate = growth_rate_liquid(sp.vol, sp.dry_vol, sp.kappa);
            sum_gn_droplet += growth_rate * sp.n;
        }
        sp.vol += growth_rate * dt;
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
double Simulation::check_valid_vol(double vol, const double dry_vol) {
    if (vol < dry_vol) {
        //std::cerr << "Found volume = " << vol/dry_vol << " * dry volume after growth." << std::endl;
        vol = dry_vol;
    }
    return vol;
}

// Coagulates the superparticles
void Simulation::coagulation() {
    const int num_sps = pop.num_sps;
    // Initialise coagulation flag vectors (std::atomic<bool> ensures thread-safe updates)
    std::vector<std::atomic<bool>> sps_coag_flag(pop.num_sps);
    for (int i = 0; i < num_sps; i++) {
        sps_coag_flag[i].store(false);
    }
    // Initialise empty vector for created superparticles
    std::vector<Superparticle> new_sps;
    // Initialise a distribution for selecting from sps vector
    std::uniform_int_distribution<int> sps_dist(0, num_sps-1);
    
    // Total possible pairings
    const double total_pairs = 0.5*num_sps*(num_sps-1);
    // Number of pairings chosen
    const int num_checks = num_sps;
    const double prob_scaling = total_pairs/num_checks;

    // Constant part of probability expression
    // Assumes all superparticles have same number density
    // Which they should for 1-to-1 coagulation
    const double const_fac = pop.n_tot/pop.num_sps * num_dt_for_coag*dt * prob_scaling;

    // Parallelise for loop
    #pragma omp parallel for schedule(dynamic) default(none) \
            shared(num_checks, sps_dist, coag_dist, const_fac, pop, sps_coag_flag, new_sps, std::cerr)
    for (int n = 0; n < num_checks; n++) {
        // Choose random pair of superparticles
        intPair pair = coag_rand_sps(sps_dist, sps_coag_flag);
        if (pair.i == -1 || pair.j == -1) {
            // Couldn't find uncoagulated superparticles; give up
            continue;
        }
        const Superparticle& sp_i = pop.sps.at(pair.i);
        const Superparticle& sp_j = pop.sps.at(pair.j);
        double beta_ij = coag_coeff(sp_i, sp_j);
        double random_number = coag_dist(global_rng());
        double P_ij = beta_ij * const_fac;
        if (P_ij >= 1) {
            #pragma omp critical
            {
                std::cerr << "Found coagulation probability P_ij = " << P_ij << ". Reduce dt for accuracy." << std::endl;
            }
        }
        if (random_number < P_ij) {
            // Coagulation occurs
            sps_coag_flag[pair.i].store(true);
            sps_coag_flag[pair.j].store(true);
            Superparticle new_sp = make_new_sp(sp_i, sp_j);
            #pragma omp critical
            {
                new_sps.push_back(new_sp);
            }
        }
    }

    // Remove coagulated superparticles
    // Create a predicate to read the flags
    auto flag_it = sps_coag_flag.begin();
    auto flag_end = sps_coag_flag.end();
    auto predicate = [&](Superparticle) {
        if (flag_it == flag_end) {
            return false;
        }
        // Returns the value at the position the iterator is pointing to
        return (*flag_it++).load();
    };

    auto new_end = std::remove_if(pop.sps.begin(), pop.sps.end(), predicate);
    pop.sps.erase(new_end, pop.sps.end());

    // Add new superparticles
    pop.sps.insert(pop.sps.end(), new_sps.begin(), new_sps.end());

    // Re-evaluate n_tot and num_sps
    pop.update_n_tot();
    pop.update_num_sps();
}

// Chooses a random pair of indices for coagulation
intPair Simulation::coag_rand_sps(std::uniform_int_distribution<>& sps_dist, std::vector<std::atomic<bool>>& sps_coag_flag) {
    intPair pair = {-1, -1};
    for (int n = 0; n < num_rand_attempts; n++) {
        int random_index = sps_dist(global_rng());
        if (!sps_coag_flag[random_index].load()) {
            pair.i = random_index;
            break;
        }
    }
    for (int n = 0; n < num_rand_attempts; n++) {
        int random_index = sps_dist(global_rng());
        if (!sps_coag_flag[random_index].load() && random_index != pair.i) {
            pair.j = random_index;
            break;
        }
    }
    return pair;
}

// Calculates the Brownian coagulation coefficient (m3 s-1) from Fuchs' interpolation formula (Seinfeld and Pandis Table 13.1)
double Simulation::coag_coeff(const Superparticle& sp_i, const Superparticle& sp_j) {
    const double vi = sp_i.vol;
    const double vj = sp_j.vol;
    const double ri = v_to_r(vi);
    const double rj = v_to_r(vj);
    double Di = diffusivity(ri);
    double Dj = diffusivity(rj);
    double ci = thermal_speed(vi, sp_i.isFrozen);
    double cj = thermal_speed(vj, sp_j.isFrozen);
    double gi = coag_g(ri, Di, ci);
    double gj = coag_g(rj, Dj, cj);
    double denom = ((ri + rj)/(ri + rj + std::sqrt(std::pow(gi, 2) + std::pow(gj, 2)))) + ((4*(Di + Dj))/(std::sqrt(std::pow(ci, 2) + std::pow(cj, 2)) * (ri + rj)));
    double beta_ij = 4*PI * (Di+Dj) * (ri + rj) / denom;
    return beta_ij;
}

// Calculates diffusivity (m2 s-1)
double Simulation::diffusivity(const double r) {
    double Cc = cscf(r);
    double D = BOLTZMANN_CONSTANT * env.get_T() * Cc / (6 * PI * env.get_air_viscosity() * r);
    return D;
}

// Calculates Cunningham slip correction factor
double Simulation::cscf(const double r) {
    double Kn = env.get_mfp_air() / r;
    double Cc = 1 + Kn * (1.257 + 0.4*std::exp(-1.1/Kn));
    return Cc;
}

// Calculates thermal speed of particles with volume v (m s-1); assumes density of pure water/ice
double Simulation::thermal_speed(const double v, const bool isFrozen) {
    double density = isFrozen ? env.get_ice_density() : env.get_water_density();
    double c = std::sqrt(8 * BOLTZMANN_CONSTANT * env.get_T() / (PI * density * v));
    return c;
}

// Calculates g for Fuchs coagulation (see Seinfeld and Pandis Table 13.1)
double Simulation::coag_g(const double r, const double D, const double c) {
    double l = 8*D/(PI*c);
    double d = 2*r;
    double g = std::sqrt(2)/(3*d*l) * (std::pow(d + l, 3) - std::pow(std::pow(d, 2) + std::pow(l, 2), 1.5)) - d;
    return g;
}

// Calculates the resulting particle properties from a coagulation and returns them in an SPTemp object
Superparticle Simulation::make_new_sp(const Superparticle& sp_i, const Superparticle& sp_j) {
    int new_ID;
    #pragma omp critical
    {
        new_ID = pop.sp_ID_count++;
    }
    double new_n = std::min(sp_i.n, sp_j.n); // should be the same
    double new_vol = sp_i.vol + sp_j.vol;
    double new_dry_vol = sp_i.dry_vol + sp_j.dry_vol;
    double new_kappa = (sp_i.vol*sp_i.kappa + sp_j.vol*sp_j.kappa)/new_vol;
    double new_ice_germs = sp_i.ice_germs + sp_j.ice_germs;
    bool isFrozen = (sp_i.isFrozen || sp_j.isFrozen);
    return Superparticle(new_ID, new_n, new_vol, new_dry_vol, new_kappa, new_ice_germs, isFrozen);
}

// Updates number of ice germs in droplets and freezes if >1
void Simulation::freezing() {
    // Update number of ice germs using Koop et al. 2000
    double u_w_i_minus_u_w_0 = 210368. + 131.438 * env.get_T() - 3.32373e6/env.get_T() - 41729.1 * std::log(env.get_T());
    double a_w_i = std::exp(u_w_i_minus_u_w_0/(IDEAL_GAS_CONSTANT * env.get_T()));
    double v_w_0 = -230.76 - 0.1478 * env.get_T() + 4099.2/env.get_T() + 48.8341 * std::log(env.get_T());
    double v_i_0 = 19.43 - 2.2e-3 * env.get_T() + 1.08e-5 * std::pow(env.get_T(), 2);
    #pragma omp parallel for
    for (Superparticle& sp : pop.sps) {
        if (!sp.isFrozen) {
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

            // Freeze if number of ice germs >= 1
            if (sp.ice_germs >= 1) {
                sp.isFrozen = true;
            }
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

    // Add distribution to n_droplet_output and n_crystal_output
    for (const Superparticle& sp : pop.sps) {
        double r = v_to_r(sp.vol);
        for (int i = 0; i < num_r_intervals_output; i++) {
            if (r >= r_output.at(i) && r < r_output.at(i+1)) {
                if (sp.isFrozen) {
                    n_crystal_output.at(i) += sp.n;
                }
                else {
                    n_droplet_output.at(i) += sp.n;
                }
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