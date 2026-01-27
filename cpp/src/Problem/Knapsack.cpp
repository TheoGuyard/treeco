#include "treeco/Problem/Knapsack.hpp"

namespace treeco {

Knapsack::Knapsack(Index numItems) : numItems_(numItems) {
  dimension_ = numItems_;
  weights_.resize(numItems_);
  for (Index i = 0; i < numItems_; ++i) {
    weights_[i] = static_cast<double>(i + 1);
  }
  capacity_ = static_cast<double>(numItems_);
}

Knapsack::Knapsack(const RealVector &weights, double capacity)
    : numItems_(weights.size()), weights_(weights), capacity_(capacity) {
  dimension_ = numItems_;
  if (weights_.size() != numItems_) {
    throw std::invalid_argument("Weights size must match number of items.");
  }
  for (const auto &w : weights_) {
    if (w <= 0.0) {
      throw std::invalid_argument("All weights must be positive.");
    }
  }
  if (capacity_ <= 0.0) {
    throw std::invalid_argument("Capacity must be positive.");
  }
}

Domain Knapsack::getCostDomain() const { return positiveOrthant(dimension_); }

std::vector<BinaryVector> Knapsack::getFeasibleSet() const {
  std::vector<BinaryVector> feasibleSet;
  Index totalSolutions = 1ULL << numItems_; // 2^numItems_
  for (Index sol = 0; sol < totalSolutions; ++sol) {
    double totalWeight = 0.0;
    BinaryVector x(numItems_, 0);
    for (Index item = 0; item < numItems_; ++item) {
      if (sol & (1ULL << item)) { // Check if item bit is set to 1
        x[item] = 1;
        totalWeight += weights_[item];
      }
    }
    if (totalWeight <= capacity_) {
      feasibleSet.push_back(x);
    }
  }
  return feasibleSet;
}

RealVector Knapsack::sampleCost(Index seed) const {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::mt19937 gen(seed);
  RealVector cost(dimension_);
  for (Index i = 0; i < dimension_; ++i) {
    cost[i] = 1.0 - dist(gen);
  }
  return cost;
}

} // namespace treeco
