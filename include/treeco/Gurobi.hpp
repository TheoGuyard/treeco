/**
 * @file Gurobi.hpp
 * @brief Gurobi solver environment management.
 * 
 * This header provides a singleton Gurobi environment that is shared
 * across all LP/MIP operations in the library.
 */

#ifndef TREECO_GUROBI_HPP
#define TREECO_GUROBI_HPP

#include <cstdlib>
#include <mutex>
#include <string>

#include "gurobi_c++.h"

/**
 * @brief Get the shared Gurobi environment.
 * 
 * Returns a singleton GRBEnv instance configured for silent operation.
 * Thread count is automatically set from SLURM_CPUS_PER_TASK if available.
 * 
 * @return Reference to the shared Gurobi environment
 */
inline GRBEnv& getGurobiEnv() {
    static std::once_flag flag;
    static GRBEnv* env = nullptr;

    const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
    int nthreads = (slurm_cpus != nullptr) ? std::stoi(slurm_cpus) : 0;

    std::call_once(flag, [nthreads] {
        env = new GRBEnv(true);
        env->set(GRB_IntParam_OutputFlag, 0);
        env->set(GRB_IntParam_Threads, nthreads);
        env->start();
    });
    
    return *env;
}

#endif // TREECO_GUROBI_HPP