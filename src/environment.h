#ifndef ENVIRONMENT
#define ENVIRONMENT

class Environment {
public:

    void initialise(std::string input_path);

    void set_env(const double current_time);

    void update_Pvap_after_growth(double delta_Pvap);

    // getters
    double get_T() {return T;}
    double get_Pvap() {return Pvap;}
    double get_T_exhaust() {return T_exhaust;}
    double get_T_ambient() {return T_ambient;}
    double get_Pvap_exhaust() {return Pvap_exhaust;}
    double get_Pvap_ambient() {return Pvap_ambient;}
    double get_P_ambient() {return P_ambient;}
    double get_r_0() {return r_0;}
    double get_eps_diffusivity() {return eps_diffusivity;}
    double get_u_0() {return u_0;}
    double get_tau_m() {return tau_m;}
    double get_Psat_l() {return Psat_l;}
    double get_Psat_i() {return Psat_i;}
    double get_S_l() {return S_l;}
    double get_S_i() {return S_i;}
    double get_air_density() {return air_density;}
    double get_air_viscosity() {return air_viscosity;}
    double get_mfp_air() {return mfp_air;}
    double get_water_density() {return water_density;}
    double get_water_density_old() {return water_density_old;}
    double get_ice_density() {return ice_density;}
    double get_ice_density_old() {return ice_density_old;}
    double get_H2O_vol_liquid() {return H2O_vol_liquid;}
    double get_H2O_vol_ice() {return H2O_vol_ice;}
    double get_water_molar_vol() {return water_molar_vol;}
    double get_ice_molar_vol() {return ice_molar_vol;}
    double get_sigma_water() {return sigma_water;}
    double get_sigma_ice() {return sigma_ice;}
    double get_vapour_thermal_speed() {return vapour_thermal_speed;}
    double get_diffusivity() {return diffusivity;}
    double get_k_air() {return k_air;}
    double get_l_v() {return l_v;}
    double get_n_sat() {return n_sat;}

private:
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

    bool first_call = true;
    
    void read_env(std::string input_path);

    double rho_w_liq(double T, double P);

    double rho_w_ice(double T, double P);
};

#endif