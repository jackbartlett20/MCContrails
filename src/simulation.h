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

    double check_valid_vol(double vol, const double dry_vol);

    void coagulation();

    double coag_coeff(const Superparticle& sp_i, const Superparticle& sp_j);

    double diffusivity(const double r);

    double cscf(const double r);

    double thermal_speed(const double v, const bool isFrozen);

    double coag_g(const double r, const double D, const double c);

    Superparticle make_new_sp(const Superparticle& sp_i, const Superparticle& sp_j);

    void freezing();

    void prepare_output();

    void output();
};

#endif