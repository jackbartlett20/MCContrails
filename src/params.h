#ifndef PARAMS
#define PARAMS

#include <string>
#include <vector>

class Species; // Forward declaration of class in population.h

class Params{
public:
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
    std::string outputDir;

    double T_exhaust;
    double T_ambient;
    double EI_vap;
    double Pvap_ambient;
    double P_ambient;
    double r_0;
    double u_0;
    double eps_diffusivity;
    double Q;
    double eta;

    std::vector<Species> species_vec;

    void read_yaml(std::string input_path);

    void read_env();

    void check_valid();

private:
    template<typename T>
    T check_and_overwrite(T var, std::string var_str, const int version=0);
};

#endif