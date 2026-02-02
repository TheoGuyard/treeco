#include "treeco/Adjacency.hpp"

namespace treeco {

Adjacency::Adjacency(const std::vector<SimplexVector>& pool)
  : pool_(pool), env_(getGurobiEnv()), model_(env_) {
  if (pool_.empty()) { throw std::invalid_argument("Pool cannot be empty"); }
  if (pool_.size() < 2) {
    throw std::invalid_argument("Pool must contain at least 2 points");
  }

  for (const auto& p : pool_) {
    if (p.size() != pool_[0].size()) {
      throw std::invalid_argument("Points must have the same dimension");
    }
  }

  build();
}

void Adjacency::build() {
  Index numPool = pool_.size();
  Index dimPool = pool_[0].size();

  // Initialize variables
  var1_.reserve(numPool);
  var2_.reserve(numPool);
  for (Index i = 0; i < numPool; ++i) {
    var1_.emplace_back(model_.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    var2_.emplace_back(
        model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
  }

  // Initialize constraints
  cstr1_.reserve(dimPool);
  cstr4_.reserve(numPool);
  for (Index j = 0; j < dimPool; ++j) {
    GRBLinExpr expr1 = 0;
    for (Index i = 0; i < numPool; ++i) { expr1 += pool_[i][j] * var2_[i]; }
    cstr1_.emplace_back(model_.addConstr(expr1 == 0.0));
  }
  GRBLinExpr expr2 = 0;
  GRBLinExpr expr3 = 0;
  for (Index i = 0; i < numPool; ++i) {
    expr2 += 1.0 * var1_[i];
    expr3 += 0.0 * var1_[i];
    cstr4_.emplace_back(model_.addConstr(var1_[i] - var2_[i] == 0.0));
  }
  cstr2_ = model_.addConstr(expr2 == 1.0);
  cstr3_ = model_.addConstr(expr3 == 1.0);

  // Initialize objective function
  GRBLinExpr obj = 0.0;
  model_.setObjective(obj, GRB_MINIMIZE);

  // Remove Gurobi outputs and set dual simplex for fast re-optimization
  model_.set("OutputFlag", "0");
  model_.set("Method", "1");
}

void Adjacency::add(Index i, Index j) {
  model_.chgCoeff(cstr2_, var1_[i], 0.0);
  model_.chgCoeff(cstr2_, var1_[j], 0.0);
  model_.chgCoeff(cstr3_, var1_[i], 1.0);
  model_.chgCoeff(cstr3_, var1_[j], 1.0);
  model_.chgCoeff(cstr4_[i], var2_[i], 1.0);
  model_.chgCoeff(cstr4_[j], var2_[j], 1.0);
}

void Adjacency::remove(Index i, Index j) {
  model_.chgCoeff(cstr2_, var1_[i], 1.0);
  model_.chgCoeff(cstr2_, var1_[j], 1.0);
  model_.chgCoeff(cstr3_, var1_[i], 0.0);
  model_.chgCoeff(cstr3_, var1_[j], 0.0);
  model_.chgCoeff(cstr4_[i], var2_[i], -1.0);
  model_.chgCoeff(cstr4_[j], var2_[j], -1.0);
}

bool Adjacency::check(Index i, Index j) {
  if (i == j) {
    throw std::invalid_argument("Index i and j must be different");
  }

  add(i, j);
  model_.optimize();
  int status = model_.get(GRB_IntAttr_Status);
  remove(i, j);

  numSolve_ += 1;

  if (status == GRB_INFEASIBLE) {
    return true;
  } else if (status == GRB_OPTIMAL) {
    return false;
  } else {
    throw std::runtime_error("Unexpected feasibility status.");
  }
}

}  // namespace treeco
