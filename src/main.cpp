#include <iostream>
#include <omp.h>
#include <chrono>
#include "simulation.h"

int main() {
    int num_threads = omp_get_num_threads();
    std::cout << "Number of threads: " << num_threads << std::endl;
    if (num_threads == 0) {
        std::cerr << "Ensure OMP_NUM_THREADS is set correctly in environment. Stopping." << std::endl;
        exit(EXIT_FAILURE);
    }
    auto start_time = std::chrono::high_resolution_clock::now();
    Simulation sim;
    sim.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Program executed successfully in " << duration.count() << " s." << std::endl;
    return 0;
}