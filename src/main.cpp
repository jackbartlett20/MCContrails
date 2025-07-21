#include <iostream>
#include <chrono>
#include "simulation.h"

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    Simulation sim;
    sim.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Program executed successfully in " << duration.count() << " s." << std::endl;
    return 0;
}