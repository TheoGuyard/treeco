#include "treeco/Problem/Explicit.hpp"

namespace treeco {

Explicit::Explicit(const std::vector<BinaryVector> &feasibleSet)
    : feasibleSet_(feasibleSet) {
  if (feasibleSet_.empty()) {
    throw std::invalid_argument("Explicit: feasible set cannot be empty");
  }
  dimension_ = feasibleSet_[0].size();
  for (const auto &x : feasibleSet_) {
    if (x.size() != dimension_) {
      throw std::invalid_argument(
          "Explicit: all vectors must have the same dimension");
    }
  }
}

std::vector<BinaryVector> Explicit::getFeasibleSet() const {
  return feasibleSet_;
}

} // namespace treeco
