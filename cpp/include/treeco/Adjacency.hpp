/**
 * @file Adjacency.hpp
 * @brief Adjacency checking for convex hull vertices.
 *
 * This header provides functionality to check whether two points in a
 * hypercube point pool are adjacent vertices of the convex hull.
 */

#ifndef TREECO_ADJACENCY_HPP
#define TREECO_ADJACENCY_HPP

#include <gurobi_c++.h>

#include <stdexcept>
#include <vector>

#include "treeco/Gurobi.hpp"
#include "treeco/Types.hpp"

namespace treeco {

/**
 * @brief Checker for adjacency between points in a hypercube point pool.
 *
 * Two points p_i and p_j from a pool P are adjacent if and only if the
 * following linear program is infeasible:
 *
 *   sum_{k != {i,j}} a_k * p_k =  a_i * p_i + a_j * p_j
 *   sum_{k != {i,j}} a_k       =  1
 *                    a_i + a_j =  1
 *                    a_k       >= 0 for all k
 *
 * This means that the segment [p_i, p_j] cannot be expressed as a convex
 * combination of other points in the pool, i.e., p_i and p_j are adjacent
 * vertices of the convex hull of P.
 *
 * For efficiency, the LP model is built once in the constructor and only
 * the coefficients are modified for each adjacency check. Dual simplex is
 * used for fast re-optimization.
 */
class Adjacency {
public:
  /**
   * @brief Construct an adjacency checker for a point pool.
   * @param pool Points in {-1,+1}^n to check adjacency from
   */
  explicit Adjacency(const std::vector<SimplexVector>& pool);

  /// Destructor
  ~Adjacency() = default;

  /**
   * @brief Check if two points in the pool are adjacent.
   * @param i Index of the first point
   * @param j Index of the second point
   * @return true if points are adjacent vertices of the convex hull
   */
  bool check(Index i, Index j);

  /// Get the number of LP solves performed
  Index getNumSolve() const { return numSolve_; }

private:
  const std::vector<SimplexVector>& pool_;  // Reference to the point pool
  Index numSolve_ = 0;                      // LP solve counter

  GRBEnv env_;                    // Gurobi environment
  GRBModel model_;                // LP model
  std::vector<GRBVar> var1_;      // Variables for coefficients
  std::vector<GRBVar> var2_;      // Variables for coefficients
  std::vector<GRBConstr> cstr1_;  // Constraint set 1
  GRBConstr cstr2_;               // Constraint 2
  GRBConstr cstr3_;               // Constraint 3
  std::vector<GRBConstr> cstr4_;  // Constraint set 4

  void build();
  void add(Index i, Index j);
  void remove(Index i, Index j);
};

}  // namespace treeco

#endif  // TREECO_ADJACENCY_HPP
