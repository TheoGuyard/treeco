/**
 * @file Adjacency.hpp
 * @brief Adjacency checker for convex hull vertices.
 */

#ifndef TREECO_ADJACENCY_HPP
#define TREECO_ADJACENCY_HPP

#include <stdexcept>
#include <vector>

#include <gurobi_c++.h>

#include "treeco/Gurobi.hpp"
#include "treeco/Types.hpp"

namespace treeco {

/**
 * @brief Checker for adjacency between vertices of the convex hull formed by 
 * a pool of points in {-1,+1}^n.
 *
 * @note Two vertices p_i and p_j of the pool are adjacent iff the system
 *
 *   sum_{k != {i,j}} a_k * p_k =  a_i * p_i + a_j * p_j
 *   sum_{k != {i,j}} a_k       =  1
 *                    a_i + a_j =  1
 *                    a_k       >= 0 for all k
 *
 * is infeasible. This means that the segment [p_i, p_j] cannot be expressed as
 * a convex combination of other points in the pool. For efficiency, the linear
 * system is built once in the constructor and only the coefficients
 * are modified for each adjacency check. Dual simplex is used for fast
 * re-optimization.
 */
class Adjacency {
public:
  
  explicit Adjacency(const std::vector<SimplexVector>& pool);

  ~Adjacency() = default;

  /**
   * @brief Check if two points in the pool are adjacent vertices.
   * @param i Index of the first point
   * @param j Index of the second point
   * @return Whether the points are adjacent vertices of the convex hull
   */
  bool check(Index i, Index j);

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
