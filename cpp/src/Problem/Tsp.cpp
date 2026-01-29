#include "treeco/Problem/Tsp.hpp"

namespace treeco {

Tsp::Tsp(Index numCities) : numCities_(numCities) {
  if (numCities_ < 3) { throw std::invalid_argument("Tsp: need at least 3 cities"); }
  dimension_ = numCities_ * (numCities_ - 1) / 2;
}

Domain Tsp::getCostDomain() const { return negativeOrthant(dimension_); }

std::vector<BinaryVector> Tsp::getFeasibleSet() const {
  std::vector<BinaryVector> tours;

  // Number of unique Hamiltonian cycles = (d-1)!/2 (undirected graph)
  Index numTours = 1;
  for (Index k = 2; k < numCities_; ++k) { numTours *= k; }
  numTours /= 2;
  tours.reserve(numTours);

  std::vector<Index> perm(numCities_ - 1);
  std::iota(perm.begin(), perm.end(), 1);
  do {
    // To avoid counting both directions of the same cycle,
    // only include tours where perm[0] < perm.back()
    // This breaks the symmetry since reversing a tour swaps these
    if (perm[0] > perm.back()) { continue; }

    BinaryVector tour(dimension_, 0);
    tour[edgeToIndex(0, perm[0])] = 1;
    for (Index k = 0; k + 1 < perm.size(); ++k) { tour[edgeToIndex(perm[k], perm[k + 1])] = 1; }
    tour[edgeToIndex(0, perm.back())] = 1;
    tours.push_back(std::move(tour));
  } while (std::next_permutation(perm.begin(), perm.end()));

  tours.shrink_to_fit();

  return tours;
}

RealVector Tsp::sampleCost(Index seed) const {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::mt19937 gen(seed);

  // Sample random city positions in the unit 2d square
  std::vector<std::pair<double, double>> positions(numCities_);
  for (Index i = 0; i < numCities_; ++i) { positions[i] = {unif(gen), unif(gen)}; }

  // Compute opposite of distances as costs
  RealVector cost(dimension_, 0.0);
  for (Index i = 0; i < numCities_; ++i) {
    for (Index j = i + 1; j < numCities_; ++j) {
      double dx = positions[i].first - positions[j].first;
      double dy = positions[i].second - positions[j].second;
      double dt = std::sqrt(dx * dx + dy * dy);
      cost[edgeToIndex(i, j)] = -dt;
    }
  }

  return cost;
}

Index Tsp::edgeToIndex(Index i, Index j) const {
  if (i >= j) { std::swap(i, j); }
  return i * numCities_ - i * (i + 1) / 2 + (j - i - 1);
}

std::pair<Index, Index> Tsp::indexToEdge(Index idx) const {
  Index i = 0;
  Index c = 0;
  while (c + (numCities_ - 1 - i) <= idx) {
    c += numCities_ - 1 - i;
    ++i;
  }
  Index j = i + 1 + (idx - c);
  return {i, j};
}

}  // namespace treeco
