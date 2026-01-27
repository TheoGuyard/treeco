/**
 * @file Tsp.hpp
 * @brief Traveling salesman problem with fixed number of cities.
 */

#ifndef TREECO_PROBLEM_TSP_HPP
#define TREECO_PROBLEM_TSP_HPP

#include <algorithm>
#include <numeric>
#include <utility>

#include "treeco/Geometry.hpp"
#include "treeco/Problem.hpp"

namespace treeco {

/**
 * @brief Traveling Salesman Problem on a complete undirected graph.
 *
 * Given d cities, the problem is to find a Hamiltonian cycle (tour visiting
 * each city exactly once) that minimizes the total distance.
 *
 * Formulation as max <c, x>:
 *   - Variables: x_ij for edges (i,j) with i < j, so d*(d-1)/2 edge variables
 *   - x_ij = 1 if edge (i,j) is in the tour
 *   - c_ij = negative distance (for minimization via maximization)
 *   - Feasible set X = all binary vectors corresponding to valid Hamiltonian
 * cycles
 *
 * Edge indexing: edge (i,j) with i < j has index i*d - i*(i+1)/2 + (j-i-1)
 */
class Tsp : public Problem {
public:
  /**
   * @brief Construct a TSP instance.
   * @param numCities Number of cities in the complete graph
   */
  explicit Tsp(Index numCities);

  /**
   * @brief Get the number of cities.
   * @return Number of cities in the graph.
   */
  Index numCities() const { return numCities_; }

  /**
   * @brief Get all feasible solutions (all hamiltonian paths).
   * @return Vector of all feasible solutions, represented as x = {x_ij} for
   * all edges (i,j) with i < j.
   */
  std::vector<BinaryVector> getFeasibleSet() const override;

  /**
   * @brief Get the cost vector domain (positive entires).
   * @return Domain of the cost vector C = {c | c_ij > 0 for all edges (i,j)}
   */
  Domain getCostDomain() const override;

  using Problem::sampleCost;

  /**
   * @brief Sample a random cost vector using a provided seed.
   * @param seed Seed for the random number generator
   * @return Cost vector for randomly places cities in [0,1)^2
   */
  RealVector sampleCost(Index seed) const override;

  /**
   * @brief Convert edge indices to variable index.
   * @param i First city (i < j)
   * @param j Second city
   * @return Variable index for edge (i,j)
   */
  Index edgeToIndex(Index i, Index j) const;

  /**
   * @brief Convert variable index to edge indices.
   * @param idx Variable index
   * @return Pair (i, j) with i < j
   */
  std::pair<Index, Index> indexToEdge(Index idx) const;

private:
  Index numCities_; // Number of cities
};

} // namespace treeco

#endif // TREECO_PROBLEM_TSP_HPP
