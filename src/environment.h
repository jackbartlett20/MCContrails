#ifndef ENVIRONMENT
#define ENVIRONMENT

class Environment {
public:
    // Variables
    double T;
    double Pvap;
    double T_exhaust;
    double T_ambient;
    double Pvap_exhaust;
    double Pvap_ambient;
    double P_ambient;
    double r_0;
    double eps_diffusivity;
    double u_0;
    double tau_m;
    double Psat_l;
    double Psat_i;
    double S_l;
    double S_i;
    double air_density;
    double air_viscosity;
    double mfp_air;
    double water_density;
    double water_density_old;
    double ice_density;
    double ice_density_old;
    double H2O_vol_liquid;
    double H2O_vol_ice;
    double water_molar_vol;
    double ice_molar_vol;
    double sigma_water;
    double sigma_ice;
    double vapour_thermal_speed;
    double diffusivity;
    double k_air;
    double l_v;
    double n_sat;

    Environment();

    void set_env(const double current_time);

private:
    // Variables
    bool first_call = true;
    
    void read_env();

    double rho_w_liq(double T, double P);

    double rho_w_ice(double T, double P);
};

#endif