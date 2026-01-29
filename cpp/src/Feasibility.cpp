#include "treeco/Feasibility.hpp"

namespace treeco {

Feasibility::Feasibility(const std::vector<TernaryVector>& pool, const Domain& fixedConstraints, double tolerance)
  : pool_(pool), tolerance_(tolerance), env_(getGurobiEnv()), model_(env_) {
  if (pool_.empty()) { throw std::invalid_argument("Pool cannot be empty"); }

  for (const auto& p : pool_) {
    if (p.size() != pool_[0].size()) { throw std::invalid_argument("Points must have the same dimension"); }
  }

  for (const auto& [w, r, b] : fixedConstraints) {
    if (w.size() != pool_[0].size()) { throw std::invalid_argument("Fixed cuts dimension mismatch."); }
  }

  if (tolerance_ < 1e-8) { throw std::invalid_argument("Tolerance must be >= 1e-8."); }

  build(fixedConstraints);
}

void Feasibility::build(const Domain& fixedConstraints) {
  Index numPool = pool_.size();
  Index dimPool = pool_[0].size();
  double slacksUb = 10.0 * tolerance_;

  // Objective expression
  GRBLinExpr objexpr_ = 0.0;

  // Variables x (dimension of the space)
  varx_.reserve(dimPool);
  for (Index j = 0; j < dimPool; ++j) {
    varx_.emplace_back(model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
  }

  // Slack variables (to handle strict inequalities)
  vars_.reserve(numPool);
  for (Index i = 0; i < numPool; ++i) { vars_.emplace_back(model_.addVar(0.0, slacksUb, 0.0, GRB_CONTINUOUS)); }

  // Add slacks to the objective
  for (Index i = 0; i < numPool; ++i) { objexpr_ += vars_[i]; }

  // Precompute linear expressions in the pool of cuts
  linexprs_.reserve(numPool);
  for (Index i = 0; i < numPool; ++i) {
    GRBLinExpr linexpr_ = 0.0;
    for (Index j = 0; j < dimPool; ++j) { linexpr_ += pool_[i][j] * varx_[j]; }
    linexpr_ += vars_[i];
    linexprs_.emplace_back(linexpr_);
  }

  // Add fixed constrains
  varsFixed_ = std::vector<GRBVar>();
  for (const auto& [a, b, r] : fixedConstraints) {
    if (r == Relation::RF) {
      throw std::runtime_error("RF relation encountered in fixed constraints.");
    } else if (r == Relation::RT) {
      continue;  // Do not add fixed constraints always satisfied
    }

    // Linear expressions for the fixed constraints
    GRBLinExpr linexprFixed = b;
    for (Index j = 0; j < dimPool; ++j) { linexprFixed += a[j] * varx_[j]; }

    // Slack for the fixed constraints if needed
    if (r == Relation::LT) {
      GRBVar varsFixed = model_.addVar(0.0, slacksUb, 0.0, GRB_CONTINUOUS);
      linexprFixed += varsFixed;
      objexpr_ += varsFixed;
      varsFixed_.push_back(varsFixed);
    } else if (r == Relation::GT) {
      GRBVar varsFixed = model_.addVar(0.0, slacksUb, 0.0, GRB_CONTINUOUS);
      linexprFixed -= varsFixed;
      objexpr_ += varsFixed;
      varsFixed_.push_back(varsFixed);
    }

    // Add the fixed constraint to the model
    if (r == Relation::LT) {
      model_.addConstr(linexprFixed, '<', 0.0);
    } else if (r == Relation::LE) {
      model_.addConstr(linexprFixed, '<', 0.0);
    } else if (r == Relation::EQ) {
      model_.addConstr(linexprFixed, '=', 0.0);
    } else if (r == Relation::GE) {
      model_.addConstr(linexprFixed, '>', 0.0);
    } else if (r == Relation::GT) {
      model_.addConstr(linexprFixed, '>', 0.0);
    } else {
      throw std::invalid_argument("Unexpected relation in fixed constraints.");
    }
  }

  // Set the objective
  model_.setObjective(objexpr_, GRB_MAXIMIZE);

  // Set Gurobi tolerance to fit tolerance_ parameter
  model_.set(GRB_DoubleParam_FeasibilityTol, 0.1 * tolerance_);

  // Remove Gurobi outputs and set dual simplex for fast re-optimization
  model_.set("OutputFlag", "0");
  model_.set("Method", "1");

  // Initialize the counts and reduced version of relations
  const std::unordered_map<Relation, Index> initialCounts = {{Relation::RT, 0}, {Relation::LT, 0}, {Relation::LE, 0},
                                                             {Relation::EQ, 0}, {Relation::GE, 0}, {Relation::GT, 0},
                                                             {Relation::RF, 0}};
  relationsCounts_.resize(numPool, initialCounts);
  relationsReduce_.resize(numPool, Relation::RT);

  // Initialize feasibility status_ as unknown
  status_ = FeasibilityStatus::UNKNOWN;
}

void Feasibility::add(const Cut& cut) {
  const auto& [i, r] = cut;
  relationsCounts_[i][r]++;
  if (relationsCounts_[i][r] == 1) {
    bool flag = reduce(i);
    if (flag) { update(i, true); }
  }
}

void Feasibility::add(const std::set<Cut>& cuts) {
  for (const auto& cut : cuts) { add(cut); }
}

void Feasibility::remove(const Cut& cut) {
  const auto& [i, r] = cut;
  if (relationsCounts_[i][r] <= 0) { throw std::runtime_error("Attempting to remove un-existing relation."); }
  relationsCounts_[i][r]--;
  if (relationsCounts_[i][r] == 0) {
    bool flag = reduce(i);
    if (flag) { update(i, false); }
  }
}

void Feasibility::remove(const std::set<Cut>& cuts) {
  for (const auto& cut : cuts) { remove(cut); }
}

bool Feasibility::reduce(Index i) {
  Relation previousRelation = relationsReduce_[i];

  std::vector<Relation> activeRelations;
  activeRelations.reserve(relationsCounts_[i].size());
  for (const auto& [r, c] : relationsCounts_[i]) {
    if (c > 0) { activeRelations.push_back(r); }
  }

  relationsReduce_[i] = reduceRelations(activeRelations);

  return (relationsReduce_[i] != previousRelation);
}

void Feasibility::update(Index i, bool isNew) {
  // New relation that must be applied to the i-th pool element of the system
  Relation r = relationsReduce_[i];

  // Warm start feasibility of the new system can be inferred when:
  // - any of the relations is infeasible, or
  // - the previous system was feasible, the new relation was already in the
  // system (hence the slack is properly set), and it is still feasible.
  // Otherwise, the feasibility of the new system is unknown.
  bool anyInfeasible =
      std::find(relationsReduce_.begin(), relationsReduce_.end(), Relation::RF) != relationsReduce_.end();
  if (anyInfeasible) {
    status_ = FeasibilityStatus::INFEASIBLE;
  } else if (status_ == FeasibilityStatus::FEASIBLE && !isNew) {
    double value = linexprs_[i].getValue();
    bool stillFeasible =
        (((r == Relation::LT) && (value <= -tolerance_)) || ((r == Relation::LE) && (value <= 0.0)) ||
         ((r == Relation::EQ) && (std::abs(value) <= tolerance_)) || ((r == Relation::GE) && (value >= 0.0)) ||
         ((r == Relation::GT) && (value >= tolerance_)) || (r == Relation::RT));
    if (stillFeasible) {
      status_ = FeasibilityStatus::FEASIBLE;
    } else {
      status_ = FeasibilityStatus::UNKNOWN;
    }
  } else {
    status_ = FeasibilityStatus::UNKNOWN;
  }

  // Update the model to match the relation
  bool alreadyInModel = indexToCstr_.find(i) != indexToCstr_.end();
  if (r == Relation::LT) {
    if (alreadyInModel) {
      indexToCstr_[i].set(GRB_CharAttr_Sense, '<');
      indexToCstr_[i].set(GRB_DoubleAttr_RHS, 0.0);
    } else {
      indexToCstr_[i] = model_.addConstr(linexprs_[i], '<', 0.0);
    }
    model_.chgCoeff(indexToCstr_[i], vars_[i], 1.0);
  } else if (r == Relation::LE) {
    if (alreadyInModel) {
      indexToCstr_[i].set(GRB_CharAttr_Sense, '<');
      indexToCstr_[i].set(GRB_DoubleAttr_RHS, 0.0);
    } else {
      indexToCstr_[i] = model_.addConstr(linexprs_[i], '<', 0.0);
    }
    model_.chgCoeff(indexToCstr_[i], vars_[i], 0.0);
  } else if (r == Relation::EQ) {
    if (alreadyInModel) {
      indexToCstr_[i].set(GRB_CharAttr_Sense, '=');
      indexToCstr_[i].set(GRB_DoubleAttr_RHS, 0.0);
    } else {
      indexToCstr_[i] = model_.addConstr(linexprs_[i], '=', 0.0);
    }
    model_.chgCoeff(indexToCstr_[i], vars_[i], 0.0);
  } else if (r == Relation::GE) {
    if (alreadyInModel) {
      indexToCstr_[i].set(GRB_CharAttr_Sense, '>');
      indexToCstr_[i].set(GRB_DoubleAttr_RHS, 0.0);
    } else {
      indexToCstr_[i] = model_.addConstr(linexprs_[i], '>', 0.0);
    }
    model_.chgCoeff(indexToCstr_[i], vars_[i], 0.0);
  } else if (r == Relation::GT) {
    if (alreadyInModel) {
      indexToCstr_[i].set(GRB_CharAttr_Sense, '>');
      indexToCstr_[i].set(GRB_DoubleAttr_RHS, 0.0);
    } else {
      indexToCstr_[i] = model_.addConstr(linexprs_[i], '>', 0.0);
    }
    model_.chgCoeff(indexToCstr_[i], vars_[i], -1.0);
  } else if (r == Relation::RT) {
    // Remove the constraint since it should always be satisfied
    if (alreadyInModel) {
      model_.remove(indexToCstr_[i]);
      indexToCstr_.erase(i);
    }
  } else if (r == Relation::RF) {
    // Remove the constraint (the model will be flagged infeasible anyway
    // due to the presence of a RF relation)
    if (alreadyInModel) {
      model_.remove(indexToCstr_[i]);
      indexToCstr_.erase(i);
    }
  } else {
    throw std::runtime_error("Unknown relation encountered.");
  }
}

bool Feasibility::check() {
  if (status_ == FeasibilityStatus::UNKNOWN) {
    model_.optimize();
    numSolve_++;

    if (model_.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
      status_ = FeasibilityStatus::INFEASIBLE;
    } else if (model_.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      status_ = FeasibilityStatus::FEASIBLE;

      // Check slacks feasibility according to tolerance_
      for (const auto& slack : vars_) {
        if (slack.get(GRB_DoubleAttr_X) < tolerance_) {
          status_ = FeasibilityStatus::INFEASIBLE;
          break;
        }
      }
      for (const auto& slack : varsFixed_) {
        if (slack.get(GRB_DoubleAttr_X) < tolerance_) {
          status_ = FeasibilityStatus::INFEASIBLE;
          break;
        }
      }
    }
  }

  if (status_ == FeasibilityStatus::FEASIBLE) {
    return true;
  } else if (status_ == FeasibilityStatus::INFEASIBLE) {
    return false;
  } else {
    throw std::runtime_error("Unexpected feasibility status.");
  }
}

}  // namespace treeco