/**
 * @file Feasibility.hpp
 * @brief Feasibility checking for linear systems of inequalities.
 *
 * This header provides functionality to incrementally check feasibility
 * of linear systems defined by cuts from a vector pool.
 */

#ifndef TREECO_FEASIBILITY_HPP
#define TREECO_FEASIBILITY_HPP

#include <gurobi_c++.h>

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "treeco/Geometry.hpp"
#include "treeco/Gurobi.hpp"

namespace treeco {

/**
 * @brief Status of feasibility for linear systems of inequalities.
 */
enum class FeasibilityStatus : int8_t {
  FEASIBLE = 0,    // System is feasible
  INFEASIBLE = 1,  // System is infeasible
  UNKNOWN = 2,     // Feasibility unknown
};

/**
 * @brief Checker for feasibility of linear systems defined by vector cuts.
 *
 * The linear system is defined dynamically by adding and removing cuts
 * from a pool of vectors V on the hypercube {-1,0,+1}^n:
 *
 *   <vi,x> <  0   for all i with Relation::LT
 *   <vi,x> <= 0   for all i with Relation::LE
 *   <vi,x> =  0   for all i with Relation::EQ
 *   <vi,x> >= 0   for all i with Relation::GE
 *   <vi,x> >  0   for all i with Relation::GT
 *
 * Strict inequalities are handled using slack variables with margin
 * maximization. Fixed constraints can be added at construction time
 * and cannot be removed later.
 */
class Feasibility {
public:
  /**
   * @brief Construct a feasibility checker.
   * @param pool Normal vectors defining possible cuts
   * @param fixedConstraints Fixed constraints included permanently
   * @param tolerance Numerical tolerance for strict inequality margins
   */
  Feasibility(const std::vector<TernaryVector>& pool, const Domain& fixedConstraints = Domain(),
              double tolerance = 1e-8);

  /// Destructor
  ~Feasibility() = default;

  /**
   * @brief Add multiple cuts to the system.
   * @param cuts Set of cuts to add
   */
  void add(const std::set<Cut>& cuts);

  /**
   * @brief Remove multiple cuts from the system.
   * @param cuts Set of cuts to remove
   */
  void remove(const std::set<Cut>& cuts);

  /**
   * @brief Add a single cut to the system.
   * @param cut The cut to add
   */
  void add(const Cut& cut);

  /**
   * @brief Remove a single cut from the system.
   * @param cut The cut to remove
   */
  void remove(const Cut& cut);

  /**
   * @brief Check feasibility of the current system.
   * @return true if feasible, false otherwise
   */
  bool check();

  /// Get the number of LP solves performed
  Index numSolve() const { return numSolve_; }

private:
  const std::vector<TernaryVector>& pool_;                 // Reference to the vector pool
  double tolerance_;                                       // Numerical tolerance
  FeasibilityStatus status_ = FeasibilityStatus::UNKNOWN;  // Current status
  Index numSolve_ = 0;                                     // LP solve counter

  GRBEnv env_;                        // Gurobi environment
  GRBModel model_;                    // LP model
  std::vector<GRBVar> varx_;          // Decision variables
  std::vector<GRBVar> vars_;          // Slack variables
  std::vector<GRBLinExpr> linexprs_;  // Linear expressions for pool vectors
  std::vector<GRBVar> varsFixed_;     // Fixed constraint variables

  std::unordered_map<Index, GRBConstr> indexToCstr_;                  // Mapping from pool index to constraint
  std::vector<std::unordered_map<Relation, Index>> relationsCounts_;  // Relation counts per pool index
  std::vector<Relation> relationsReduce_;                             // Reduced relations per pool index

  void build(const Domain& fixedConstraints);
  bool reduce(Index i);
  void update(Index i, bool isNew);
};

}  // namespace treeco

#endif  // TREECO_FEASIBILITY_HPP