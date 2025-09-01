#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <calculisto/iapws/r6_inverse.hpp>
#include <calculisto/iapws/r10.hpp>
#include "environment.h"
#include "constants.h"

Environment::Environment() {
    read_env();
    T = T_exhaust;
    Pvap = Pvap_exhaust;
    double x_m = r_0 * std::sqrt(2/eps_diffusivity);
    tau_m = x_m / u_0;
    set_env(0);
}

// Reads input file with variables relevant to Environment
void Environment::read_env() {
    std::ifstream file("input/environment.in");
    if (!file.is_open()) {
        std::cerr << "Error opening environment.in" << std::endl;
        exit(EXIT_FAILURE);
    }
    file >> T_exhaust;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> T_ambient;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> Pvap_exhaust;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> Pvap_ambient;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> P_ambient;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> r_0;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> u_0;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> eps_diffusivity;
    file.close();

    // Check valid
    if (T_exhaust <= 0) {
        std::cerr << "Error: Read in exhaust temperature of " << T_exhaust << "." << std::endl;
        std::cerr << "Exhaust temperature must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (T_ambient <= 0) {
        std::cerr << "Error: Read in ambient temperature of " << T_ambient << "." << std::endl;
        std::cerr << "Ambient temperature must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (Pvap_exhaust < 0) {
        std::cerr << "Error: Read in exhaust vapour pressure of " << Pvap_exhaust << "." << std::endl;
        std::cerr << "Exhaust vapour pressure must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (Pvap_ambient < 0) {
        std::cerr << "Error: Read in ambient vapour pressure of " << Pvap_ambient << "." << std::endl;
        std::cerr << "Ambient vapour pressure must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (P_ambient <= 0) {
        std::cerr << "Error: Read in ambient pressure of " << P_ambient << "." << std::endl;
        std::cerr << "Ambient pressure must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_0 <= 0) {
        std::cerr << "Error: Read in radius of jet exhaust of " << r_0 << "." << std::endl;
        std::cerr << "Radius of jet exhaust must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (u_0 <= 0) {
        std::cerr << "Error: Read in aircraft speed of " << u_0 << "." << std::endl;
        std::cerr << "Aircraft speed must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (eps_diffusivity <= 0) {
        std::cerr << "Error: Read in turbulent diffusivity of " << eps_diffusivity << "." << std::endl;
        std::cerr << "Turbulent diffusivity must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Sets environmental variables according to current time
void Environment::set_env(const double current_time) {
    // Adiabatic plume mixing
    double dilution_factor;
    if (current_time <= tau_m) {
        dilution_factor = 1;
    }
    else {
        dilution_factor = std::pow(tau_m/current_time, 0.9);
    }

    double mixing_grad = (Pvap - Pvap_ambient)/(T - T_ambient);
    double T_new = T_ambient + (T_exhaust - T_ambient) * dilution_factor;
    // Pvap can't be calculated the same way because it is also updated by growth
    // However, it does decrease like temperature due to mixing
    Pvap += mixing_grad * (T_new - T);
    T = T_new;

    // Buck equations
    Psat_l = 6.1121e2*std::exp((18.678 - (T-273.15)/234.5) * ((T-273.15)/(T-16.01)));
    Psat_i = 6.1115e2*std::exp((23.036 - (T-273.15)/333.7) * ((T-273.15)/(T+6.67)));

    // Saturation ratios
    S_l = Pvap/Psat_l;
    S_i = Pvap/Psat_i;

    // Density of air (kg m-3)
    air_density = P_ambient * AIR_MOLAR_MASS / (IDEAL_GAS_CONSTANT * T);

    // Viscosity of air (N s m-2) - Sutherland equation
    air_viscosity = 1.716e-5 * std::pow(T/273, 1.5) * 384/(T+111);

    // Mean free path of air (m) - assumes effective collision diameter 0.37 nm
    mfp_air = (BOLTZMANN_CONSTANT*T)/(std::sqrt(2)*PI*std::pow(0.37e-9, 2)*P_ambient);

    // Water and ice density (kg m-3)
    if (first_call) {
        water_density_old = rho_w_liq(T, P_ambient);
        ice_density_old = rho_w_ice(T, P_ambient);
    }
    else {
        water_density_old = water_density;
        ice_density_old = ice_density;
    }
    first_call = false;
    water_density = rho_w_liq(T, P_ambient);
    ice_density = rho_w_ice(T, P_ambient);

    // H2O molecular vol in liquid (m3)
    H2O_vol_liquid = H2O_MOLECULAR_MASS / water_density;

    // H2O molecular vol in ice (m3)
    H2O_vol_ice = H2O_MOLECULAR_MASS / ice_density;

    // H2O molar volume in liquid (m3 mol-1)
    water_molar_vol = H2O_vol_liquid * AVOGADRO_CONSTANT;

    // H2O molar volume in ice (m3 mol-1)
    ice_molar_vol = H2O_vol_ice * AVOGADRO_CONSTANT;

    // Water surface tension (N m-1) - IAPWS R1-76(2014)
    sigma_water = 235.8e-3 * std::pow(1-T/647.096, 1.256) * (1 - 0.625*(1-T/647.096));

    // Ice surface tension (N m-1) - Pruppacher and Klett
    sigma_ice = 0.1;

    // Thermal speed of water vapour (m s-1)
    vapour_thermal_speed = std::sqrt(8 * BOLTZMANN_CONSTANT * T / (PI*H2O_MOLECULAR_MASS));

    // Diffusivity (m2 s-1) - Pruppacher and Klett eq. 13-3
    diffusivity = 2.11e-5 * std::pow(T/273.15, 1.94) * (101325 / P_ambient);

    // Thermal conductivity of air (J m-1 s-1 K-1) - Seinfeld and Pandis eq. 17.54
    k_air = 1e-3 * (4.39 + 0.071 * T);

    // Specific latent heat of vaporisation of water (J kg-1)
    l_v = 2.501e6 - 2.37e3 * (T - 273.15);

    // H2O number concentration (m-3) - for crystal growth; ideally not needed here
    n_sat = AVOGADRO_CONSTANT * Pvap / (IDEAL_GAS_CONSTANT * T);
}

// Calculates density of liquid water (kg m-3) using IAPWS R6-95(2018)
double Environment::rho_w_liq(double T, double P) {
    double rho = calculisto::iapws::r6_inverse::density_pt(P, std::max(T, 235.));
    return rho;
}

// Calculates density of ice (kg m-3) using IAPWS R10-06(2009)
double Environment::rho_w_ice(double T, double P) {
    double rho = calculisto::iapws::r10::density_pt(P, T);
    return rho;
}