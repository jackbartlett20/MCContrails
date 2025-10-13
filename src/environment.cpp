#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <calculisto/iapws/r6_inverse.hpp>
#include <calculisto/iapws/r10.hpp>
#include "environment.h"
#include "constants.h"
#include "params.h"

void Environment::initialise(Params& params) {
    init_vars(params);
    T = T_exhaust;
    Pvap = Pvap_exhaust;
    double x_m = r_0 * std::sqrt(2/eps_diffusivity);
    tau_m = x_m / u_0;
    set_env(0);
}

// Copies relevant variables from Params object to self
void Environment::init_vars(Params& params) {
    T_exhaust       = params.T_exhaust;
    T_ambient       = params.T_ambient;
    Pvap_exhaust    = params.Pvap_exhaust;
    Pvap_ambient    = params.Pvap_ambient;
    P_ambient       = params.P_ambient;
    r_0             = params.r_0;
    u_0             = params.u_0;
    eps_diffusivity = params.eps_diffusivity;
}

// Sets environmental variables according to current time
void Environment::set_env(const double current_time) {
    // Adiabatic plume mixing from Kärcher et al. (2015)
    if (first_call) {
        dilution_factor_old = (current_time <= tau_m) ? 1 : std::pow(tau_m/current_time, 0.9);
    }
    else {
        dilution_factor_old = dilution_factor;
    }
    dilution_factor = (current_time <= tau_m) ? 1 : std::pow(tau_m/current_time, 0.9);

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
    if (first_call) {
        air_density_old = P_ambient * AIR_MOLAR_MASS / (IDEAL_GAS_CONSTANT * T);
    }
    else {
        air_density_old = air_density;
    }
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

    first_call = false;
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

// Updates vapour pressure in Environment after growth
void Environment::update_Pvap_after_growth(double delta_Pvap) {
    Pvap += delta_Pvap;
}