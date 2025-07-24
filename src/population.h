#ifndef POPULATION
#define POPULATION

#include <vector>

class Species {
public:
    // Variables
    double n;
    double GMR;
    double GSD;
    double f_dry;
    double kappa;

    Species(double n, double r_mean, double SD, double f_dry, double kappa);
};

class Superparticle {
public:
    // Variables
    int ID;
    double n;
    double vol;
    double f_dry;
    double kappa;
    double ice_germs;
    bool isFrozen;

    Superparticle(int ID, double n, double vol, double f_dry, double kappa, double ice_germs, bool isFrozen);
};

// Used in coagulation to temporarily store the doubles for new superparticles
class SPTemp {
public:
    // Variables
    double n;
    double vol;
    double f_dry;
    double kappa;
    double ice_germs;

    SPTemp(double n, double vol, double f_dry, double kappa, double ice_germs);
};

class Population {
public:
    // Variables
    int num_sps;
    double n_tot;
    int max_sps;
    int num_r_choices;
    double f_dry_tol;
    int sp_ID_count;
    std::vector<Superparticle> droplet_sps;
    std::vector<Superparticle> crystal_sps;

    void assign();
    void update_n_tot();
    void update_num_sps();

private:
    void read_population();
    std::vector<Species> read_species();
    std::vector<double> choose_vs(int sps_for_species, std::vector<double> logr_range, std::vector<double> weights);
};

#endif