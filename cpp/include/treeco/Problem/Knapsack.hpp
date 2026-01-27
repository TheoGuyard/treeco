/**
 * @file Knapsack.hpp
 * @brief Knapsack Problem with fixed weights and capacities.
 */

#ifndef TREECO_PROBLEM_KNAPSACK_HPP
#define TREECO_PROBLEM_KNAPSACK_HPP

#include <algorithm>
#include <numeric>
#include <utility>

#include "treeco/Geometry.hpp"
#include "treeco/Problem.hpp"

namespace treeco {

/**
 * @brief Knapsack Problem with fixed weights and capacities.
 *
 * Given n items, the problem is to select a subset of items to maximize the
 * total value without exceeding the capacity.
 *
 * Formulation as max <c, x> s.t.  <w, x> <= W where:
 *   - Variables: x_i for items i, binary variables indicating item selection
 *   - c_i = value of item i
 *   - w_i = weight of item i
 *   - W = total knapsack capacity
 */
class Knapsack : public Problem {
public:
  /**
   * @brief Construct a Knapsack instance with weights {1, 2, ..., n} and
   * capacity n.
   * @param numItems Number of items in the knapsack
   */
  explicit Knapsack(Index numItems);

  /**
   * @brief Construct a Knapsack instance.
   * @param weights Weights of the items
   * @param capacity Total capacity of the knapsack
   */
  explicit Knapsack(const RealVector &weights, double capacity);

  /**
   * @brief Get the number of items.
   * @return Number of items in the knapsack.
   */
  Index numItems() const { return numItems_; }

  /**
   * @brief Get the weights of the items.
   * @return Vector of item weights.
   */
  RealVector weights() const { return weights_; }

  /**
   * @brief Get the knapsack capacity.
   * @return Capacity of the knapsack.
   */
  double capacity() const { return capacity_; }

  /**
   * @brief Get all feasible solutions (all valid item selections).
   * @return Vector of all feasible solutions, represented as x = {x_i} for
   * all items i.
   */
  std::vector<BinaryVector> getFeasibleSet() const override;

  /**
   * @brief Get the cost vector domain (positive entires).
   * @return Domain of the cost vector C = {c | c_i > 0 for all items i}
   */
  Domain getCostDomain() const override;

  using Problem::sampleCost;

  /**
   * @brief Sample a random cost vector using a provided seed.
   * @param seed Seed for the random number generator
   * @return Cost vector with entries uniformly distributed in (0, 1]
   */
  RealVector sampleCost(Index seed) const override;

private:
  Index numItems_;     // Number of items
  RealVector weights_; // Weights of the items
  double capacity_;    // Total capacity of the knapsack
};

} // namespace treeco

#endif // TREECO_PROBLEM_KNAPSACK_HPP
