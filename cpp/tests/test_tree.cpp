#include <gtest/gtest.h>

#include <sstream>
#include <treeco/Dynprog.hpp>
#include <treeco/Geometry.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Tree.hpp>
#include <treeco/Types.hpp>
#include <treeco/Voronoi.hpp>

using namespace treeco;

class TreeTest : public ::testing::Test {
protected:
  std::vector<SimplexVector> points;
  std::unique_ptr<Voronoi> voronoi;

  void SetUp() override {
    // Use TSP with 4 cities as test case
    Tsp tsp(4);
    points = scaleBinarySet(tsp.getFeasibleSet());
    voronoi = std::make_unique<Voronoi>(points);
    voronoi->build();
  }
};

// ============================================================================
// Construction Tests
// ============================================================================

TEST_F(TreeTest, Construction) {
  Tree tree;
  EXPECT_FALSE(tree.isBuilt());
}

TEST_F(TreeTest, ConstructionDynprogNotRun) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  EXPECT_THROW(tree.synthetize(dynprog), std::runtime_error);
}

// ============================================================================
// Synthetize Tests
// ============================================================================

TEST_F(TreeTest, SynthetizeBasic) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();

  tree.synthetize(dynprog);

  EXPECT_TRUE(tree.isBuilt());
  EXPECT_GT(tree.size(), 0);
}

TEST_F(TreeTest, SynthetizeWithParams) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();

  std::ostringstream oss;
  TreeParams params;
  params.verbose = false;
  params.outputStream = &oss;

  tree.synthetize(dynprog, params);

  EXPECT_TRUE(tree.isBuilt());
}

TEST_F(TreeTest, SynthetizeStats) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();

  tree.synthetize(dynprog);

  const TreeStats& stats = tree.stats();
  EXPECT_TRUE(stats.isBuilt);
  EXPECT_GE(stats.buildTime, 0.0);
}

// ============================================================================
// Clear Tests
// ============================================================================

TEST_F(TreeTest, Clear) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  EXPECT_TRUE(tree.isBuilt());

  tree.clear();

  EXPECT_FALSE(tree.isBuilt());
  EXPECT_EQ(tree.size(), 0);
}

// ============================================================================
// Tree Structure Tests
// ============================================================================

TEST_F(TreeTest, NodesAccessor) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  const std::vector<Node>& nodes = tree.nodes();
  EXPECT_FALSE(nodes.empty());
  EXPECT_EQ(nodes.size(), tree.size());
}

TEST_F(TreeTest, NodeAccessor) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  // Access root node
  const Node& root = tree.node(0);
  EXPECT_EQ(root.depth, 0);
  EXPECT_NE(root.type, NodeType::UNDEFINED);
}

TEST_F(TreeTest, TreeDimensions) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  EXPECT_GT(tree.size(), 0);
  EXPECT_GE(tree.width(), 1);
  EXPECT_GE(tree.depth(), 0);
}

// ============================================================================
// Node Type Tests
// ============================================================================

TEST_F(TreeTest, LeafNodesHavePoints) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  for (const Node& node : tree.nodes()) {
    if (node.type == NodeType::LEAF) {
      EXPECT_FALSE(node.pointsIds.empty());
      EXPECT_EQ(node.splitId, INVALID_INDEX);
      EXPECT_TRUE(node.children.empty());
    }
  }
}

TEST_F(TreeTest, InternalNodesHaveSplits) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  for (const Node& node : tree.nodes()) {
    if (node.type == NodeType::NODE) {
      EXPECT_TRUE(node.pointsIds.empty());
      EXPECT_NE(node.splitId, INVALID_INDEX);
      EXPECT_FALSE(node.children.empty());
    }
  }
}

// ============================================================================
// Stream Operator Tests
// ============================================================================

TEST_F(TreeTest, StreamOperator) {
  Tree tree;
  Dynprog dynprog(*voronoi);
  dynprog.run();
  tree.synthetize(dynprog);

  std::ostringstream oss;
  oss << tree;
  EXPECT_FALSE(oss.str().empty());
}
