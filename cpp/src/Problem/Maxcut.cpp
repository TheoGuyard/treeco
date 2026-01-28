#include "treeco/Problem/Maxcut.hpp"

namespace treeco {

Maxcut::Maxcut(Index numVertices) : numVertices_(numVertices) {
  if (numVertices_ < 2) {
    throw std::invalid_argument("Maxcut: need at least 2 vertices");
  }
  dimension_ = numVertices_ * (numVertices_ - 1) / 2;
}

Domain Maxcut::getCostDomain() const { return positiveOrthant(dimension_); }

std::vector<BinaryVector> Maxcut::getFeasibleSet() const {

  std::vector<BinaryVector> cuts;
  cuts.reserve(1ULL << (numVertices_ - 1));

  Index numPartitions = 1ULL << (numVertices_ - 1);
  for (Index mask = 0; mask < numPartitions; ++mask) {
    BinaryVector partition(numVertices_, 0);
    partition[0] = 0; // Fix vertex 0 to partition 0
    for (Index v = 1; v < numVertices_; ++v) {
      partition[v] = (mask >> (v - 1)) & 1;
    }
    BinaryVector cut(dimension_, 0);
    for (Index i = 0; i < numVertices_; ++i) {
      for (Index j = i + 1; j < numVertices_; ++j) {
        if (partition[i] != partition[j]) {
          cut[edgeToIndex(i, j)] = 1;
        }
      }
    }
    cuts.push_back(cut);
  }

  cuts.shrink_to_fit();

  return cuts;
}

RealVector Maxcut::sampleCost(Index seed) const {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::mt19937 gen(seed);

  // Sample random vertices positions in the unit 2d square
  std::vector<std::pair<double, double>> positions(numVertices_);
  for (Index i = 0; i < numVertices_; ++i) {
    positions[i] = {unif(gen), unif(gen)};
  }

  // Compute distances as costs
  RealVector cost(dimension_, 0.0);
  for (Index i = 0; i < numVertices_; ++i) {
    for (Index j = i + 1; j < numVertices_; ++j) {
      double dx = positions[i].first - positions[j].first;
      double dy = positions[i].second - positions[j].second;
      double dt = std::sqrt(dx * dx + dy * dy);
      cost[edgeToIndex(i, j)] = dt;
    }
  }

  return cost;
}

Index Maxcut::edgeToIndex(Index i, Index j) const {
  if (i >= j) {
    std::swap(i, j);
  }
  return i * numVertices_ - i * (i + 1) / 2 + (j - i - 1);
}

std::pair<Index, Index> Maxcut::indexToEdge(Index idx) const {
  Index i = 0;
  Index c = 0;
  while (c + (numVertices_ - 1 - i) <= idx) {
    c += numVertices_ - 1 - i;
    ++i;
  }
  Index j = i + 1 + (idx - c);
  return {i, j};
}

} // namespace treeco
