#include <iostream>
#include <string>
#include <cstdlib>
#include <type_traits>
#include <yaml-cpp/yaml.h>
#include "params.h"
#include "population.h"

// Reads yaml input file
void Params::read_yaml(std::string input_path) {
    YAML::Node input_file = YAML::LoadFile(input_path);
    
    // Simulation
    int_time               = input_file["simulation"]["int_time"].as<double>();
    dt                     = input_file["simulation"]["dt"].as<double>();
    num_writes             = input_file["simulation"]["num_writes"].as<int>();
    r_output_min           = input_file["simulation"]["r_output_min"].as<double>();
    r_output_max           = input_file["simulation"]["r_output_max"].as<double>();
    num_r_intervals_output = input_file["simulation"]["num_r_intervals_output"].as<int>();
    max_sps                = input_file["simulation"]["max_sps"].as<int>();
    num_r_choices          = input_file["simulation"]["num_r_choices"].as<int>();
    min_S_l                = input_file["simulation"]["min_S_l"].as<double>();
    do_coagulation         = input_file["simulation"]["do_coagulation"].as<int>();
    num_dt_for_coag        = input_file["simulation"]["num_dt_for_coag"].as<int>();
    rng_seed_read          = input_file["simulation"]["rng_seed"].as<unsigned long long>();
    outputDir              = input_file["simulation"]["MCC_output_dir"].as<std::string>();

    // Environment
    T_exhaust       = input_file["environment"]["T_exhaust"].as<double>();
    T_ambient       = input_file["environment"]["T_ambient"].as<double>();
    EI_vap          = input_file["environment"]["EI_vap"].as<double>();
    Pvap_ambient    = input_file["environment"]["Pvap_ambient"].as<double>();
    P_ambient       = input_file["environment"]["P_ambient"].as<double>();
    r_0             = input_file["environment"]["r_0"].as<double>();
    u_0             = input_file["environment"]["u_0"].as<double>();
    eps_diffusivity = input_file["environment"]["eps_diffusivity"].as<double>();
    Q               = input_file["environment"]["Q"].as<double>();
    eta             = input_file["environment"]["eta"].as<double>();

    // Iterate over species
    for (const auto& speciesNode : input_file["exhaust species"]) {
        double EI    = speciesNode["EI"].as<double>();
        double GMR   = speciesNode["GMR"].as<double>();
        double GSD   = speciesNode["GSD"].as<double>();
        double f_dry = speciesNode["f_dry"].as<double>();
        double kappa = speciesNode["kappa"].as<double>();
        species_vec.push_back(Species(EI, GMR, GSD, f_dry, kappa));
    }
}

// Looks for variables set in the environment and overwrites with their values
void Params::read_env() {
    // Simulation
    int_time               = check_and_overwrite<double>(int_time, "int_time");
    dt                     = check_and_overwrite<double>(dt, "dt");
    num_writes             = check_and_overwrite<int>(num_writes, "num_writes");
    r_output_min           = check_and_overwrite<double>(r_output_min, "r_output_min");
    r_output_max           = check_and_overwrite<double>(r_output_max, "r_output_max");
    num_r_intervals_output = check_and_overwrite<int>(num_r_intervals_output, "num_r_intervals_output");
    max_sps                = check_and_overwrite<int>(max_sps, "max_sps");
    num_r_choices          = check_and_overwrite<int>(num_r_choices, "num_r_choices");
    min_S_l                = check_and_overwrite<double>(min_S_l, "min_S_l");
    do_coagulation         = check_and_overwrite<int>(do_coagulation, "do_coagulation");
    num_dt_for_coag        = check_and_overwrite<int>(num_dt_for_coag, "num_dt_for_coag");
    rng_seed_read          = check_and_overwrite<unsigned long long>(rng_seed_read, "rng_seed");
    outputDir              = check_and_overwrite<std::string>(outputDir, "MCC_output_dir");

    // Environment
    T_exhaust       = check_and_overwrite<double>(T_exhaust, "T_exhaust");
    T_ambient       = check_and_overwrite<double>(T_ambient, "T_ambient");
    EI_vap          = check_and_overwrite<double>(EI_vap, "EI_vap");
    Pvap_ambient    = check_and_overwrite<double>(Pvap_ambient, "Pvap_ambient");
    P_ambient       = check_and_overwrite<double>(P_ambient, "P_ambient");
    r_0             = check_and_overwrite<double>(r_0, "r_0");
    u_0             = check_and_overwrite<double>(u_0, "u_0");
    eps_diffusivity = check_and_overwrite<double>(eps_diffusivity, "eps_diffusivity");
    Q               = check_and_overwrite<double>(Q, "Q");
    eta             = check_and_overwrite<double>(eta, "eta");

    // Iterate over species, uses index i to determine which version to look for
    for (std::vector<Species>::size_type i = 0; i < species_vec.size(); i++) {
        Species& species = species_vec.at(i);
        int version = i+1;
        species.EI     = check_and_overwrite<double>(species.EI, "n", version);
        species.GMR   = check_and_overwrite<double>(species.GMR, "GMR", version);
        species.GSD   = check_and_overwrite<double>(species.GSD, "GSD", version);
        species.f_dry = check_and_overwrite<double>(species.f_dry, "f_dry", version);
        species.kappa = check_and_overwrite<double>(species.kappa, "kappa", version);
    }
}

/*
Looks for specific variable in environment and overwrites if available.
var is the variable, var_str is the string to look for in the environment,
and version is a number to append to the species variables e.g. n1 instead of n.
For other variables, version=0.
*/
template<typename T>
T Params::check_and_overwrite(T var, std::string var_str, const int version) {
    if (version != 0) {
        var_str.append(std::to_string(version));
    }
    char* c_str = std::getenv(var_str.c_str()); // C style string
    if (c_str != nullptr) {
        // Variable exists in environment
        std::string cpp_str(c_str); // C++ style string
        if constexpr (std::is_same_v<T, std::string>) {
            var = cpp_str;
        }
        else if constexpr (std::is_same_v<T, int>) {
            var = std::stoi(cpp_str);
        }
        else if constexpr (std::is_same_v<T, float>) {
            var = std::stof(cpp_str);
        }
        else if constexpr (std::is_same_v<T, double>) {
            var = std::stod(cpp_str);
        }
        else if constexpr (std::is_same_v<T, unsigned long long>) {
            var = std::stoull(cpp_str);
        }
        std::cout << "Variable read from environment: " << var_str << "=" << var << std::endl;
    }
    return var;
}

// Checks parameters are valid
void Params::check_valid() {
    // Simulation
    if (int_time <= 0) {
        std::cerr << "Error: Read in integration time of " << int_time << "." << std::endl;
        std::cerr << "Integration time must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (dt <= 0) {
        std::cerr << "Error: Read in time step of " << dt << "." << std::endl;
        std::cerr << "Time step must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_writes < 0) {
        std::cerr << "Error: Read in " << num_writes << " writes." << std::endl;
        std::cerr << "Number of writes must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_output_min <= 0) {
        std::cerr << "Error: Read in smallest radius in output spectrum of " << r_output_min << "." << std::endl;
        std::cerr << "Smallest radius in output spectrum must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_output_max <= r_output_min) {
        std::cerr << "Error: Read in largest radius in output spectrum of " << r_output_max << "." << std::endl;
        std::cerr << "Largest radius in output spectrum must be greater than smallest radius in output spectrum. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_r_intervals_output <= 0) {
        std::cerr << "Error: Read in " << num_r_intervals_output << " intervals in output spectrum." << std::endl;
        std::cerr << "Number of intervals in output spectrum must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (max_sps <= 0) {
        std::cerr << "Error: Read in total number of superparticles of " << max_sps << "." << std::endl;
        std::cerr << "Total number of superparticles must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_r_choices <= 0) {
        std::cerr << "Error: Read in " << num_r_choices << " radii to choose from during initialisation." << std::endl;
        std::cerr << "Number of radii to choose from must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (min_S_l < 0) {
        std::cerr << "Error: Read in minimum saturation ratio of " << min_S_l << "." << std::endl;
        std::cerr << "Minimum saturation ratio must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (do_coagulation != 0 && do_coagulation != 1) {
        std::cerr << "Error: Read in \"do coagulation?\" of" << do_coagulation << "." << std::endl;
        std::cerr << "Enter 1 for coagulation, 0 for no coagulation. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (num_dt_for_coag < 0) {
        std::cerr << "Error: Read in number of time steps for coagulation of " << num_dt_for_coag << "." << std::endl;
        std::cerr << "Number of time steps for coagulation must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Environment
    if (T_exhaust <= 0) {
        std::cerr << "Error: Read in exhaust temperature of " << T_exhaust << "." << std::endl;
        std::cerr << "Exhaust temperature must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (T_ambient <= 0) {
        std::cerr << "Error: Read in ambient temperature of " << T_ambient << "." << std::endl;
        std::cerr << "Ambient temperature must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (EI_vap < 0) {
        std::cerr << "Error: Read in water vapour emissions index of " << EI_vap << "." << std::endl;
        std::cerr << "Water vapour emissions index must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (Pvap_ambient < 0) {
        std::cerr << "Error: Read in ambient vapour pressure of " << Pvap_ambient << "." << std::endl;
        std::cerr << "Ambient vapour pressure must be >= 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (P_ambient <= 0) {
        std::cerr << "Error: Read in ambient pressure of " << P_ambient << "." << std::endl;
        std::cerr << "Ambient pressure must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (r_0 <= 0) {
        std::cerr << "Error: Read in radius of jet exhaust of " << r_0 << "." << std::endl;
        std::cerr << "Radius of jet exhaust must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (u_0 <= 0) {
        std::cerr << "Error: Read in aircraft speed of " << u_0 << "." << std::endl;
        std::cerr << "Aircraft speed must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (eps_diffusivity <= 0) {
        std::cerr << "Error: Read in turbulent diffusivity of " << eps_diffusivity << "." << std::endl;
        std::cerr << "Turbulent diffusivity must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (Q <= 0) {
        std::cerr << "Error: Read in specific combustion heat of " << Q << "." << std::endl;
        std::cerr << "Specific combustion heat must be > 0. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (eta <= 0 || eta > 1) {
        std::cerr << "Error: Read in propulsion efficiency of " << eta << "." << std::endl;
        std::cerr << "Propulsion efficiency must be between 0 and 1. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Species
    for (Species& species : species_vec) {
        if (species.EI < 0) {
            std::cerr << "Error: Read in emissions index of " << species.EI << "." << std::endl;
            std::cerr << "Emissions index should be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (species.GMR <= 0) {
            std::cerr << "Error: Read in geometric mean radius of " << species.GMR << "." << std::endl;
            std::cerr << "Geometric mean radius must be > 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (species.f_dry < 0 || species.f_dry > 1) {
            std::cerr << "Error: Read in dry fraction of " << species.f_dry << "." << std::endl;
            std::cerr << "Dry fraction must be between 0 and 1. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (species.kappa < 0) {
            std::cerr << "Error: Read in hygroscopicity of " << species.kappa << "." << std::endl;
            std::cerr << "Hygroscopicity must be >= 0. Stopping." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}