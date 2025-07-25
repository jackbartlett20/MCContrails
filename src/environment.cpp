#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <cmath>
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
}

// Sets environmental variables according to current time
void Environment::set_env(const double current_time) {
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
        water_density_old = rho_w_liq(T);
        ice_density_old = rho_w_ice(T);
    }
    else {
        water_density_old = water_density;
        ice_density_old = ice_density;
    }
    first_call = false;
    water_density = rho_w_liq(T);
    ice_density = rho_w_ice(T);

    // H2O molecular vol in liquid (m3)
    H2O_vol_liquid = H2O_MOLECULAR_MASS / water_density;

    // H2O molecular vol in ice (m3)
    H2O_vol_ice = H2O_MOLECULAR_MASS / ice_density;

    // Water surface tension (N m-1) - IAPWS
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

    // H2O number concentration (m-3) - REMOVE WHEN CRYSTAL GROWTH UPDATED
    n_sat = AVOGADRO_CONSTANT * Pvap / (IDEAL_GAS_CONSTANT * T);
}

// Calculates density of liquid water (kg m-3) - Sippola and Taskinen (2018)
double Environment::rho_w_liq(double T) {
    double rho_0 = 1007.853; // kg m-3
    double T_c = 228; // K
    double A = 3.9744e-4;
    double B = 1.6785e-3;
    double C = -7.816510e-4;
    double eps;
    if (T > T_c) {
        eps = T/T_c - 1;
    }
    else {
        eps = 0;
    }
    double rho = rho_0 * std::exp(-T_c * (A + B*eps + 2*C*std::pow(eps, 0.5)));
    return rho;
}

// Calculates approximate density of ice (kg m-3)
double Environment::rho_w_ice(double T) {
    double rho = 917 - 0.13*(std::min(T-273.15, 0.));
    return rho;
}