#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>
#include "simulation.h"

int main(int argc, char* argv[]) {
    // Check input file
    std::string input_path;
    if (argc > 1) {
        // Input file path provided in the command line
        input_path = argv[1];
    }
    else {
        // No argument; use default path to input.
        input_path = "../input.yaml";
    }
    std::ifstream file(input_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open input file at '" << input_path << "'. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    file.close();
    std::cout << "Using input file: " << input_path << std::endl;

    // Check thread count
    int num_threads = 0;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    std::cout << "Number of threads: " << num_threads << std::endl;
    if (num_threads == 0) {
        std::cerr << "Ensure OMP_NUM_THREADS is set correctly in environment. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Do simulation
    auto start_time = std::chrono::high_resolution_clock::now();
    Simulation sim;
    sim.run(input_path);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Program executed successfully in " << duration.count() << " s." << std::endl;
    return 0;
}