#include <gtest/gtest.h>
#include <treeco/Feasibility.hpp>

using namespace treeco;

// ============================================================================
// Test fixture with common setup
// ============================================================================

class FeasibilityTest : public ::testing::Test {
protected:
  // Simple 2D pool: directions including axis-aligned (with 0s)
  std::vector<TernaryVector> pool2D = {
      {1, 0},  // x1
      {0, 1},  // x2
      {1, 1},  // x1 + x2
      {1, -1}, // x1 - x2
  };

  // 3D pool
  std::vector<TernaryVector> pool3D = {
      {1, 0, 0}, // x1
      {0, 1, 0}, // x2
      {0, 0, 1}, // x3
      {1, 1, 0}, // x1 + x2
      {1, 0, 1}, // x1 + x3
      {0, 1, 1}, // x2 + x3
      {1, 1, 1}, // x1 + x2 + x3
  };
};

// ============================================================================
// Constructor tests
// ============================================================================

TEST_F(FeasibilityTest, ConstructorValid) {
  EXPECT_NO_THROW((Feasibility(pool2D)));
}

TEST_F(FeasibilityTest, ConstructorEmptyPool) {
  std::vector<TernaryVector> emptyPool;
  EXPECT_THROW((Feasibility(emptyPool)), std::invalid_argument);
}

TEST_F(FeasibilityTest, ConstructorDimensionMismatch) {
  std::vector<TernaryVector> badPool = {
      {1, 0}, {1, 0, 1} // Different dimension
  };
  EXPECT_THROW((Feasibility(badPool)), std::invalid_argument);
}

TEST_F(FeasibilityTest, ConstructorWithFixedConstraints) {
  Domain fixed = {
      {{1.0, 0.0}, 0.0, Relation::GE} // x1 >= 0
  };
  EXPECT_NO_THROW((Feasibility(pool2D, fixed)));
}

TEST_F(FeasibilityTest, ConstructorFixedConstraintDimensionMismatch) {
  Domain fixed = {
      {{1.0, 0.0, 0.0}, 0.0, Relation::GE} // Wrong dimension
  };
  EXPECT_THROW((Feasibility(pool2D, fixed)), std::invalid_argument);
}

TEST_F(FeasibilityTest, ConstructorToleranceTooSmall) {
  EXPECT_THROW((Feasibility(pool2D, {}, 1e-10)), std::invalid_argument);
}

// ============================================================================
// Basic feasibility tests
// ============================================================================

TEST_F(FeasibilityTest, EmptySystemFeasible) {
  Feasibility feas(pool2D);
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, SingleConstraintFeasible) {
  Feasibility feas(pool2D);

  // x1 >= 0 is feasible
  feas.add(Cut(0, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, TwoCompatibleConstraints) {
  Feasibility feas(pool2D);

  // x1 >= 0 and x2 >= 0 (first quadrant)
  feas.add(Cut(0, Relation::GE));
  feas.add(Cut(1, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, ContradictoryConstraintsInfeasible) {
  Feasibility feas(pool2D);

  // x1 >= 0 and x1 <= 0 reduces to x1 = 0, still feasible
  feas.add(Cut(0, Relation::GE));
  feas.add(Cut(0, Relation::LE));
  EXPECT_TRUE(feas.check()); // x1 = 0 is valid
}

TEST_F(FeasibilityTest, StrictContradictionInfeasible) {
  Feasibility feas(pool2D);

  // x1 > 0 and x1 < 0 is infeasible
  feas.add(Cut(0, Relation::GT));
  feas.add(Cut(0, Relation::LT));
  EXPECT_FALSE(feas.check());
}

TEST_F(FeasibilityTest, StrictAndNonStrictContradiction) {
  Feasibility feas(pool2D);

  // x1 > 0 and x1 <= 0 is infeasible
  feas.add(Cut(0, Relation::GT));
  feas.add(Cut(0, Relation::LE));
  EXPECT_FALSE(feas.check());
}

// ============================================================================
// Add/Remove operations
// ============================================================================

TEST_F(FeasibilityTest, AddRemoveSingleCut) {
  Feasibility feas(pool2D);

  Cut cut(0, Relation::GT);
  feas.add(cut);
  EXPECT_TRUE(feas.check());

  feas.remove(cut);
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, AddRemoveMultipleCuts) {
  Feasibility feas(pool2D);

  std::set<Cut> cuts = {Cut(0, Relation::GE), Cut(1, Relation::GE)};

  feas.add(cuts);
  EXPECT_TRUE(feas.check());

  feas.remove(cuts);
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, RemoveRestoresFeasibility) {
  Feasibility feas(pool2D);

  // Add contradictory constraints
  feas.add(Cut(0, Relation::GT));
  feas.add(Cut(0, Relation::LT));
  EXPECT_FALSE(feas.check());

  // Remove one constraint, should become feasible
  feas.remove(Cut(0, Relation::LT));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, RemoveNonExistentThrows) {
  Feasibility feas(pool2D);

  // Removing a cut that was never added should throw
  EXPECT_THROW(feas.remove(Cut(0, Relation::GE)), std::runtime_error);
}

TEST_F(FeasibilityTest, AddSameCutMultipleTimes) {
  Feasibility feas(pool2D);

  Cut cut(0, Relation::GE);
  feas.add(cut);
  feas.add(cut); // Add same cut twice
  EXPECT_TRUE(feas.check());

  // Need to remove twice
  feas.remove(cut);
  EXPECT_TRUE(feas.check()); // Still one remaining
  feas.remove(cut);
  EXPECT_TRUE(feas.check()); // Now fully removed
}

// ============================================================================
// Relation reduction tests
// ============================================================================

TEST_F(FeasibilityTest, RelationReductionLEAndGE) {
  Feasibility feas(pool2D);

  // x1 <= 0 and x1 >= 0 should reduce to x1 = 0
  feas.add(Cut(0, Relation::LE));
  feas.add(Cut(0, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, RelationReductionLTAndLE) {
  Feasibility feas(pool2D);

  // x1 < 0 and x1 <= 0 should reduce to x1 < 0
  feas.add(Cut(0, Relation::LT));
  feas.add(Cut(0, Relation::LE));
  EXPECT_TRUE(feas.check());
}

// ============================================================================
// Complex systems
// ============================================================================

TEST_F(FeasibilityTest, BoundedRegion2D) {
  Feasibility feas(pool2D);

  // Define a bounded region: x1 >= 0, x2 >= 0, x1 + x2 <= 0
  // This means x1 >= 0, x2 >= 0, x1 + x2 <= 0
  // Only solution is x1 = x2 = 0
  feas.add(Cut(0, Relation::GE)); // x1 >= 0
  feas.add(Cut(1, Relation::GE)); // x2 >= 0
  feas.add(Cut(2, Relation::LE)); // x1 + x2 <= 0
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, InfeasibleTriangle) {
  Feasibility feas(pool2D);

  // x1 > 0, x2 > 0, x1 + x2 < 0 is infeasible
  feas.add(Cut(0, Relation::GT)); // x1 > 0
  feas.add(Cut(1, Relation::GT)); // x2 > 0
  feas.add(Cut(2, Relation::LT)); // x1 + x2 < 0
  EXPECT_FALSE(feas.check());
}

TEST_F(FeasibilityTest, Feasibility3D) {
  Feasibility feas(pool3D);

  // First octant: x1 >= 0, x2 >= 0, x3 >= 0
  feas.add(Cut(0, Relation::GE));
  feas.add(Cut(1, Relation::GE));
  feas.add(Cut(2, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, Feasibility3DWithDiagonal) {
  Feasibility feas(pool3D);

  // x1 >= 0, x2 >= 0, x3 >= 0, x1 + x2 + x3 <= 0
  // Only x = 0 is valid
  feas.add(Cut(0, Relation::GE));
  feas.add(Cut(1, Relation::GE));
  feas.add(Cut(2, Relation::GE));
  feas.add(Cut(6, Relation::LE)); // x1 + x2 + x3 <= 0
  EXPECT_TRUE(feas.check());
}

// ============================================================================
// Fixed constraints tests
// ============================================================================

TEST_F(FeasibilityTest, FixedConstraintRestrictsSpace) {
  // Fixed constraint: x1 >= 1
  Domain fixed = {
      {{1.0, 0.0}, -1.0, Relation::GE} // x1 - 1 >= 0, i.e., x1 >= 1
  };
  Feasibility feas(pool2D, fixed);

  // Adding x1 <= 0 should make it infeasible
  feas.add(Cut(0, Relation::LE)); // x1 <= 0
  EXPECT_FALSE(feas.check());
}

TEST_F(FeasibilityTest, FixedConstraintCompatible) {
  // Fixed constraint: x1 >= 0
  Domain fixed = {{{1.0, 0.0}, 0.0, Relation::GE}};
  Feasibility feas(pool2D, fixed);

  // Adding x1 >= 0 from pool should still be feasible
  feas.add(Cut(0, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, FixedConstraintWithEquality) {
  // Fixed constraint: x1 = 0
  Domain fixed = {{{1.0, 0.0}, 0.0, Relation::EQ}};
  Feasibility feas(pool2D, fixed);

  EXPECT_TRUE(feas.check());

  // Adding x1 > 0 should be infeasible
  feas.add(Cut(0, Relation::GT));
  EXPECT_FALSE(feas.check());
}

// ============================================================================
// Edge cases
// ============================================================================

TEST_F(FeasibilityTest, EqualityConstraint) {
  Feasibility feas(pool2D);

  // x1 = 0 is feasible
  feas.add(Cut(0, Relation::EQ));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, MultipleEqualitiesOnSameIndex) {
  Feasibility feas(pool2D);

  // Adding EQ multiple times
  feas.add(Cut(0, Relation::EQ));
  feas.add(Cut(0, Relation::EQ));
  EXPECT_TRUE(feas.check());

  feas.remove(Cut(0, Relation::EQ));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, SequentialAddRemove) {
  Feasibility feas(pool2D);

  // Build up constraints then tear down
  feas.add(Cut(0, Relation::GE));
  EXPECT_TRUE(feas.check());

  feas.add(Cut(1, Relation::GE));
  EXPECT_TRUE(feas.check());

  feas.add(Cut(0, Relation::LT)); // Now x1 >= 0 and x1 < 0 -> infeasible
  EXPECT_FALSE(feas.check());

  feas.remove(Cut(0, Relation::LT));
  EXPECT_TRUE(feas.check());

  feas.remove(Cut(1, Relation::GE));
  EXPECT_TRUE(feas.check());

  feas.remove(Cut(0, Relation::GE));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, RTRelationAlwaysTrue) {
  Feasibility feas(pool2D);

  // RT (always true) should not affect feasibility
  feas.add(Cut(0, Relation::RT));
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, RFRelationAlwaysFalse) {
  Feasibility feas(pool2D);

  // RF (always false) should make system infeasible
  feas.add(Cut(0, Relation::RF));
  EXPECT_FALSE(feas.check());

  // Removing should restore feasibility
  feas.remove(Cut(0, Relation::RF));
  EXPECT_TRUE(feas.check());
}

// ============================================================================
// Stress/Performance tests
// ============================================================================

TEST_F(FeasibilityTest, ManyAddsAndRemoves) {
  Feasibility feas(pool2D);

  // Add all GE constraints
  for (Index i = 0; i < pool2D.size(); ++i) {
    feas.add(Cut(i, Relation::GE));
  }
  EXPECT_TRUE(feas.check());

  // Remove all
  for (Index i = 0; i < pool2D.size(); ++i) {
    feas.remove(Cut(i, Relation::GE));
  }
  EXPECT_TRUE(feas.check());
}

TEST_F(FeasibilityTest, AlternatingFeasibility) {
  Feasibility feas(pool2D);

  for (int iter = 0; iter < 5; ++iter) {
    // Make infeasible
    feas.add(Cut(0, Relation::GT));
    feas.add(Cut(0, Relation::LT));
    EXPECT_FALSE(feas.check());

    // Make feasible again
    feas.remove(Cut(0, Relation::GT));
    feas.remove(Cut(0, Relation::LT));
    EXPECT_TRUE(feas.check());
  }
}
