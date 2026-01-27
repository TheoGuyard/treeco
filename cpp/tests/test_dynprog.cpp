#include <gtest/gtest.h>
#include <treeco/Dynprog.hpp>
#include <treeco/Geometry.hpp>
#include <treeco/Problem/Maxcut.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Types.hpp>
#include <treeco/Voronoi.hpp>

using namespace treeco;

class DynprogTest : public ::testing::Test {
protected:
  std::vector<SimplexVector> smallPoints;
  std::vector<SimplexVector> mediumPoints;

  void SetUp() override {
    // Small test case: 4 points (simple)
    smallPoints = {{1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}};

    // Medium test case: TSP with n=4
    Tsp tsp(4);
    mediumPoints = scaleBinarySet(tsp.getFeasibleSet());
  }
};

TEST_F(DynprogTest, BasicConstruction) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // Should not throw
  Dynprog dp(voronoi);
  EXPECT_EQ(dp.status(), DynprogStatus::INVALID);
}

TEST_F(DynprogTest, SmallCaseBinary) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.branching = Branching::BINARY;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // With 4 faces, binary branching needs at least ceil(log2(4)) = 2 levels
  EXPECT_GE(depth, 2);
  EXPECT_LT(depth, MAX_DEPTH);
}

TEST_F(DynprogTest, SmallCaseTernary) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.branching = Branching::TERNARY;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // With 4 faces, ternary branching needs at least ceil(log3(4)) = 2 levels
  EXPECT_GE(depth, 2);
  EXPECT_LT(depth, MAX_DEPTH);
}

TEST_F(DynprogTest, MediumCaseBinary) {
  Voronoi voronoi(mediumPoints);
  VoronoiParams vparams;
  vparams.verbose = false;
  voronoi.build(vparams);

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.branching = Branching::BINARY;
  params.timeLimit = 10.0;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // Should find some valid depth
  EXPECT_LT(depth, MAX_DEPTH);

  // TSP(4) has 3 feasible tours, so at least ceil(log2(3)) = 2 levels
  EXPECT_GE(depth, 2);
}

TEST_F(DynprogTest, IterativeVsNonIterative) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // Exhaustive (single pass with all splits)
  Dynprog dp1(voronoi);
  DynprogParams params1;
  params1.verbose = false;
  params1.exploration = Exploration::EXHAUSTIVE;
  dp1.run(params1);
  Index depth1 = dp1.stats().optimalDepth;

  // Iterative
  Dynprog dp2(voronoi);
  DynprogParams params2;
  params2.verbose = false;
  params2.exploration = Exploration::ITERATIVE;
  dp2.run(params2);
  Index depth2 = dp2.stats().optimalDepth;

  // Both should give the same optimal depth
  EXPECT_EQ(depth1, depth2);
}

TEST_F(DynprogTest, DifferentScoringStrategies) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  std::vector<SplitScoring> strategies = {
      SplitScoring::VARIANCE, SplitScoring::ENTROPY, SplitScoring::MINMAX,
      SplitScoring::NONE, SplitScoring::RANDOM};

  Index baseDepth = MAX_DEPTH;
  for (auto strategy : strategies) {
    Dynprog dp(voronoi);
    DynprogParams params;
    params.verbose = false;
    params.splitScoring = strategy;
    dp.run(params);
    Index depth = dp.stats().optimalDepth;

    if (baseDepth == MAX_DEPTH) {
      baseDepth = depth;
    }
    // All strategies should find the same optimal depth
    EXPECT_EQ(depth, baseDepth);
  }
}

TEST_F(DynprogTest, DifferentLowerBoundStrategies) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // FIXED strategy
  Dynprog dp1(voronoi);
  DynprogParams params1;
  params1.verbose = false;
  params1.lowerBounding = LowerBounding::FIXED;
  dp1.run(params1);
  Index depth1 = dp1.stats().optimalDepth;

  // BACKTRACK strategy
  Dynprog dp2(voronoi);
  DynprogParams params2;
  params2.verbose = false;
  params2.lowerBounding = LowerBounding::BACKTRACK;
  dp2.run(params2);
  Index depth2 = dp2.stats().optimalDepth;

  // Both should find the same optimal depth
  EXPECT_EQ(depth1, depth2);
}

TEST_F(DynprogTest, PrecomputeVsOnFly) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // On-the-fly
  Dynprog dp1(voronoi);
  DynprogParams params1;
  params1.verbose = false;
  params1.positioning = Positioning::ONLINE;
  dp1.run(params1);
  Index depth1 = dp1.stats().optimalDepth;

  // Precompute
  Dynprog dp2(voronoi);
  DynprogParams params2;
  params2.verbose = false;
  params2.positioning = Positioning::PRECOMPUTE;
  dp2.run(params2);
  Index depth2 = dp2.stats().optimalDepth;

  // Both should find the same optimal depth
  EXPECT_EQ(depth1, depth2);
}

TEST_F(DynprogTest, WithInputConstraints) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // Add constraint: c[0] >= 0
  Domain constraints;
  constraints.emplace_back(RealVector{-1.0, 0.0, 0.0}, 0.0, Relation::LE);

  Dynprog dp(voronoi, constraints);
  DynprogParams params;
  params.verbose = false;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // Should still find a valid depth
  EXPECT_LT(depth, MAX_DEPTH);
}

TEST_F(DynprogTest, TimeLimit) {
  Voronoi voronoi(mediumPoints);
  VoronoiParams vparams;
  vparams.verbose = false;
  voronoi.build(vparams);

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.timeLimit = 0.001; // Very short time limit

  // Should return without hanging
  dp.run(params);

  // Just verify it completed - may or may not have found optimal
  EXPECT_TRUE(true);
}

TEST_F(DynprogTest, StatisticsPopulated) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;

  dp.run(params);

  const DynprogStats &stats = dp.stats();
  EXPECT_GT(stats.numStates, 0);
  EXPECT_GT(stats.numStatesClosed, 0);
  EXPECT_GE(stats.runTime, 0.0);
  EXPECT_LT(stats.optimalDepth, MAX_DEPTH);
}

TEST_F(DynprogTest, SingleFaceCase) {
  // Single point case - trivial Voronoi
  std::vector<SimplexVector> singlePoint = {{1, 1}};

  // This will fail in Voronoi constructor (needs at least 2 points)
  // So we test with 2 points where one dominates
  std::vector<SimplexVector> twoPoints = {{1, 1}, {-1, -1}};

  Voronoi voronoi(twoPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // With 2 faces, need at least 1 split
  EXPECT_GE(depth, 1);
}

TEST_F(DynprogTest, StatusAfterRun) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  EXPECT_EQ(dp.status(), DynprogStatus::INVALID);

  DynprogParams params;
  params.verbose = false;
  dp.run(params);

  // After successful run, status should be OPTIMAL or SUBOPTIMAL
  EXPECT_NE(dp.status(), DynprogStatus::INVALID);
}

TEST_F(DynprogTest, GreedyExploration) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.exploration = Exploration::GREEDY;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // Greedy may not find optimal, but should find something valid
  EXPECT_LT(depth, MAX_DEPTH);
}

TEST_F(DynprogTest, SplitSelectionSampling) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  params.splitSelection = SplitSelection::SAMPLING;
  params.randomSeed = 42;

  dp.run(params);
  Index depth = dp.stats().optimalDepth;

  // Should find a valid depth
  EXPECT_LT(depth, MAX_DEPTH);
}

TEST_F(DynprogTest, RootIdAccessor) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  Dynprog dp(voronoi);
  DynprogParams params;
  params.verbose = false;
  dp.run(params);

  Index rootId = dp.rootId();
  EXPECT_NE(rootId, INVALID_INDEX);

  // Should be able to access root state
  const State &rootState = dp.state(rootId);
  EXPECT_EQ(rootState.depth(), 0);
}

TEST_F(DynprogTest, FilterChecksOption) {
  Voronoi voronoi(smallPoints);
  voronoi.build();

  // With filter checks
  Dynprog dp1(voronoi);
  DynprogParams params1;
  params1.verbose = false;
  params1.filterChecks = true;
  dp1.run(params1);
  Index depth1 = dp1.stats().optimalDepth;

  // Without filter checks
  Dynprog dp2(voronoi);
  DynprogParams params2;
  params2.verbose = false;
  params2.filterChecks = false;
  dp2.run(params2);
  Index depth2 = dp2.stats().optimalDepth;

  // Both should find the same optimal depth
  EXPECT_EQ(depth1, depth2);
}
