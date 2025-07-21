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
    double mfp_air;
    double sigma_water;
    double vapour_thermal_speed;
    double diffusivity;
    double k_air;
    double l_v;
    double n_sat;

    Environment();

    void set_env(const double current_time);

private:
    void read_env();
};

#endif