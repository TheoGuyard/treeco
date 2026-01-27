#include "treeco/Problem.hpp"

namespace treeco {

RealVector Problem::sampleCost(Index seed) const {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  RealVector cost(dimension_);
  std::mt19937 gen(seed);
  for (Index i = 0; i < dimension_; ++i) {
    cost[i] = dist(gen);
  }
  return cost;
}

RealVector Problem::sampleCost() const {
  static thread_local std::mt19937 gen{std::random_device{}()};
  Index seed = gen();
  return sampleCost(seed);
}

} // namespace treeco
