#ifndef POPULATION
#define POPULATION

#include <vector>

class Params; // Forward declaration of class in params.h

class Species {
public:
    double EI;
    double GMR;
    double GSD;
    double f_dry;
    double kappa;

    Species(double EI, double r_mean, double SD, double f_dry, double kappa);
};

class Superparticle {
public:
    int ID;
    double n;
    double vol;
    double dry_vol;
    double kappa;
    double ice_germs;
    bool isFrozen;

    Superparticle(int ID, double n, double vol, double dry_vol, double kappa, double ice_germs, bool isFrozen);
};

class Population {
public:
    int num_sps;
    double n_tot;
    int max_sps;
    int num_r_choices;
    int sp_ID_count;
    std::vector<Superparticle> sps;

    void assign(Params& params, const double initial_n_to_EI_ratio);
    void update_n_tot();
    void update_num_sps();

private:
    std::vector<double> choose_vs(int sps_for_species, std::vector<double> logr_range, std::vector<double> weights);
};

#endif