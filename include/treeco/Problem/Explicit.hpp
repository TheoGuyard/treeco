/**
 * @file Explicit.hpp
 * @brief Problem with explicitly provided feasible set.
 */

#ifndef TREECO_PROBLEM_EXPLICIT_HPP
#define TREECO_PROBLEM_EXPLICIT_HPP

#include <vector>

#include "treeco/Problem.hpp"

namespace treeco {

/**
 * @brief Optimization problem with an explicitly provided feasible set.
 * 
 * Use this class when you have a pre-computed list of all feasible solutions.
 */
class Explicit : public Problem {
public:
    /**
     * @brief Construct from an explicit feasible set.
     * @param feasibleSet Vector of all feasible binary solutions
     */
    explicit Explicit(const std::vector<BinaryVector>& feasibleSet);

    /// Get all feasible solutions
    std::vector<BinaryVector> getFeasibleSet() const override;

private:
    std::vector<BinaryVector> feasibleSet_; // The feasible solutions
};

} // namespace treeco

#endif // TREECO_PROBLEM_EXPLICIT_HPP
