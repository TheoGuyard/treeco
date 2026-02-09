#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <treeco/Problem.hpp>
#include <treeco/Problem/Explicit.hpp>
#include <treeco/Problem/Maxcut.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Types.hpp>
#include <unordered_set>

using namespace treeco;

// ============================================================================
// Explicit Tests
// ============================================================================

TEST(ExplicitProblemTest, BasicConstruction) {
  std::vector<BinaryVector> feasibleSet = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 1}};

  Explicit problem(feasibleSet);

  EXPECT_EQ(problem.dimension(), 3);
  EXPECT_EQ(problem.getFeasibleSet().size(), 4);
}

TEST(ExplicitProblemTest, EmptySetThrows) {
  std::vector<BinaryVector> emptySet;
  EXPECT_THROW((Explicit(emptySet)), std::invalid_argument);
}

TEST(ExplicitProblemTest, InconsistentDimensionsThrow) {
  std::vector<BinaryVector> badSet = {
      {0, 0, 0}, {1, 0}  // Different dimension
  };
  EXPECT_THROW((Explicit(badSet)), std::invalid_argument);
}

TEST(ExplicitProblemTest, SampleCost) {
  std::vector<BinaryVector> feasibleSet = {{0, 0}, {1, 1}};
  Explicit problem(feasibleSet);

  RealVector cost = problem.sampleCost();

  EXPECT_EQ(cost.size(), 2);
  for (double c : cost) {
    EXPECT_GE(c, 0.0);
    EXPECT_LE(c, 1.0);
  }
}

// ============================================================================
// Tsp Tests
// ============================================================================

TEST(TSPTest, BasicConstruction) {
  Tsp tsp(4);

  // 4 cities -> 4*3/2 = 6 edges
  EXPECT_EQ(tsp.numCities(), 4);
  EXPECT_EQ(tsp.dimension(), 6);
}

TEST(TSPTest, TooFewCitiesThrows) { EXPECT_THROW(Tsp(2), std::invalid_argument); }

TEST(TSPTest, EdgeIndexing) {
  Tsp tsp(4);

  EXPECT_EQ(tsp.edgeToIndex(0, 1), 0);
  EXPECT_EQ(tsp.edgeToIndex(0, 2), 1);
  EXPECT_EQ(tsp.edgeToIndex(0, 3), 2);
  EXPECT_EQ(tsp.edgeToIndex(1, 2), 3);
  EXPECT_EQ(tsp.edgeToIndex(1, 3), 4);
  EXPECT_EQ(tsp.edgeToIndex(2, 3), 5);

  EXPECT_EQ(tsp.indexToEdge(0), std::make_pair(Index(0), Index(1)));
  EXPECT_EQ(tsp.indexToEdge(3), std::make_pair(Index(1), Index(2)));
  EXPECT_EQ(tsp.indexToEdge(5), std::make_pair(Index(2), Index(3)));
}

TEST(TSPTest, FeasibleSetSize) {
  Tsp tsp3(3);
  EXPECT_EQ(tsp3.getFeasibleSet().size(), 1);

  Tsp tsp4(4);
  EXPECT_EQ(tsp4.getFeasibleSet().size(), 3);

  Tsp tsp5(5);
  EXPECT_EQ(tsp5.getFeasibleSet().size(), 12);
}

TEST(TSPTest, ToursAreValid) {
  Tsp tsp(4);
  auto tours = tsp.getFeasibleSet();

  for (const auto& tour : tours) {
    int numEdges = std::accumulate(tour.begin(), tour.end(), 0);
    EXPECT_EQ(numEdges, 4);

    std::vector<int> degree(4, 0);
    for (Index idx = 0; idx < tour.size(); ++idx) {
      if (tour[idx] == 1) {
        auto [i, j] = tsp.indexToEdge(idx);
        degree[i]++;
        degree[j]++;
      }
    }
    for (int d : degree) { EXPECT_EQ(d, 2); }
  }
}

TEST(TSPTest, ToursAreUnique) {
  Tsp tsp(4);
  auto tours = tsp.getFeasibleSet();

  std::set<BinaryVector> uniqueTours(tours.begin(), tours.end());
  EXPECT_EQ(uniqueTours.size(), tours.size());
}

// ============================================================================
// Maxcut Tests
// ============================================================================

TEST(MaxCutTest, BasicConstruction) {
  Maxcut maxcut(4);
  EXPECT_EQ(maxcut.numVertices(), 4);
  EXPECT_EQ(maxcut.dimension(), 6);
}

TEST(MaxCutTest, TooFewVerticesThrows) { EXPECT_THROW(Maxcut(1), std::invalid_argument); }

TEST(MaxCutTest, EdgeIndexing) {
  Maxcut maxcut(4);

  EXPECT_EQ(maxcut.edgeToIndex(0, 1), 0);
  EXPECT_EQ(maxcut.edgeToIndex(0, 2), 1);
  EXPECT_EQ(maxcut.edgeToIndex(2, 3), 5);

  EXPECT_EQ(maxcut.indexToEdge(0), std::make_pair(Index(0), Index(1)));
  EXPECT_EQ(maxcut.indexToEdge(5), std::make_pair(Index(2), Index(3)));
}

TEST(MaxCutTest, FeasibleSetSize) {
  Maxcut maxcut2(2);
  EXPECT_EQ(maxcut2.getFeasibleSet().size(), 2);

  Maxcut maxcut3(3);
  EXPECT_EQ(maxcut3.getFeasibleSet().size(), 4);

  Maxcut maxcut4(4);
  EXPECT_EQ(maxcut4.getFeasibleSet().size(), 8);
}

TEST(MaxCutTest, CutsAreValid) {
  Maxcut maxcut(4);
  auto cuts = maxcut.getFeasibleSet();

  BinaryVector emptyCut(6, 0);
  EXPECT_NE(std::find(cuts.begin(), cuts.end(), emptyCut), cuts.end());

  BinaryVector starCut = {1, 1, 1, 0, 0, 0};
  EXPECT_NE(std::find(cuts.begin(), cuts.end(), starCut), cuts.end());
}

TEST(MaxCutTest, CutsAreUnique) {
  Maxcut maxcut(4);
  auto cuts = maxcut.getFeasibleSet();

  std::set<BinaryVector> uniqueCuts(cuts.begin(), cuts.end());
  EXPECT_EQ(uniqueCuts.size(), cuts.size());
}

TEST(MaxCutTest, SampleCostInRange) {
  Maxcut maxcut(5);

  RealVector cost = maxcut.sampleCost();

  EXPECT_EQ(cost.size(), maxcut.dimension());
  for (double c : cost) {
    EXPECT_GE(c, 0.0);
    EXPECT_LE(c, std::sqrt(2.0));
  }
}
