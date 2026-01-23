/**
 * @file Problem.hpp
 * @brief Abstract base class for optimization problems.
 * 
 * This header defines the Problem interface for linear optimization problems
 * of the form max <c, x> subject to x in X.
 */

#ifndef TREECO_PROBLEM_PROBLEM_HPP
#define TREECO_PROBLEM_PROBLEM_HPP

#include <random>
#include <vector>

#include "treeco/Types.hpp"

namespace treeco {

/**
 * @brief Abstract base class for optimization problems.
 * 
 * Represents a linear optimization problem of the form:
 *   max <c, x>  subject to  x in X
 * 
 * where c is a real-valued cost vector and X is a set of binary vectors.
 * Subclasses implement specific problem types by defining the feasible set X.
 */
class Problem {
public:
    virtual ~Problem() = default;

    /// Get the dimension of the problem (length of vectors x and c)
    Index dimension() const { return dimension_; }

    /**
     * @brief Get all feasible solutions.
     * @return Vector of all binary vectors in the feasible set X
     */
    virtual std::vector<BinaryVector> getFeasibleSet() const = 0;

    /**
     * @brief Get the cost vector domain.
     * @return Domain of the cost vector as linear constraints.
     * @note When not overloaded, returns a domain corresponding to the whole
     * ambient space (no constraints).
     */
    Domain getCostDomain() const;

    /**
     * @brief Sample a random cost vector using a provided seed.
     * @param seed Seed for the random number generator
     * @return Cost vector with entries uniformly distributed in [0, 1)
     */
    RealVector sampleCost(Index seed) const;

    /**
     * @brief Sample a random cost vector.
     * @return Cost vector with entries uniformly distributed in [0, 1)
     */
    RealVector sampleCost() const;

protected:
    Index dimension_ = INVALID_INDEX; // Problem dimension
};

} // namespace treeco

#endif // TREECO_PROBLEM_PROBLEM_HPP
