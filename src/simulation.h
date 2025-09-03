#ifndef SIMULATION
#define SIMULATION

#include <vector>
#include "environment.h"
#include "population.h"

class Simulation {
public:

    void run(std::string input_path);

private:
    double int_time;
    double dt;
    int num_writes;
    double r_output_min;
    double r_output_max;
    int num_r_intervals_output;
    int max_sps;
    int num_r_choices;
    double min_S_l;
    int do_coagulation;
    int num_dt_for_coag;
    unsigned long long rng_seed_read;
    double write_interval;
    double next_write_time;
    double coag_interval;
    double next_coag_time;
    double current_time;
    Environment env;
    Population pop;
    std::vector<double> r_output;
    std::vector<double> r_m_output;
    std::vector<double> n_droplet_output;
    std::vector<double> dndlogr_droplet_output;
    std::vector<double> n_crystal_output;
    std::vector<double> dndlogr_crystal_output;
    double dlogr_output;
    bool first_write;

    void read_simulation(std::string input_path);

    void update_water_vol();

    void growth();

    double growth_rate_liquid(const double v, const double v_dry, const double kappa);

    double growth_rate_crystal(const double v);

    double check_valid_vol(double vol, double dry_vol);

    void coagulation();

    double coag_coeff(double vi, double vj);

    double diffusivity(double r);

    double cscf(double r);

    double thermal_speed(double v);

    double coag_g(double r, double D, double c);

    SPTemp new_props(Superparticle sp_i, Superparticle sp_j);

    void freezing();

    void prepare_output();

    void output();
};

#endif