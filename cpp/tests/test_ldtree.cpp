#include <gtest/gtest.h>
#include <treeco/Geometry.hpp>
#include <treeco/LDTree.hpp>
#include <treeco/Problem/Explicit.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Types.hpp>

#include <cstdio>
#include <fstream>
#include <sstream>

using namespace treeco;

class LDTreeTest : public ::testing::Test {
protected:
  std::vector<BinaryVector> points;
  std::string tempPointsFile;
  std::string tempOutputFile;

  void SetUp() override {
    // Use TSP with 4 cities as test case
    Tsp tsp(4);
    points = tsp.getFeasibleSet();

    // Create temporary file for file-based tests
    tempPointsFile = "/tmp/treeco_test_ldtree_points_" +
                     std::to_string(std::rand()) + ".txt";
    tempOutputFile =
        "/tmp/treeco_test_ldtree_output_" + std::to_string(std::rand()) + ".c";

    // Write points to file
    std::ofstream outfile(tempPointsFile);
    outfile << "points " << points[0].size() << " " << points.size() << "\n";
    for (const auto &point : points) {
      for (size_t i = 0; i < point.size(); ++i) {
        outfile << static_cast<int>(point[i]);
        if (i < point.size() - 1)
          outfile << " ";
      }
      outfile << "\n";
    }
    outfile.close();
  }

  void TearDown() override {
    std::remove(tempPointsFile.c_str());
    std::remove(tempOutputFile.c_str());
  }
};

// ============================================================================
// Construction Tests
// ============================================================================

TEST_F(LDTreeTest, ConstructFromPoints) {
  LDTree ldtree(points);

  EXPECT_EQ(ldtree.voronoi().numPoints(), points.size());
}

TEST_F(LDTreeTest, ConstructFromFile) {
  LDTree ldtree(tempPointsFile);

  EXPECT_EQ(ldtree.voronoi().numPoints(), points.size());
}

TEST_F(LDTreeTest, ConstructWithEmptyDomain) {
  Domain emptyDomain;
  LDTree ldtree(points, emptyDomain);

  EXPECT_TRUE(ldtree.domain().empty());
}

TEST_F(LDTreeTest, ConstructFromFileNotFound) {
  EXPECT_THROW(LDTree("/nonexistent/path.txt"), std::runtime_error);
}

// ============================================================================
// Build Tests
// ============================================================================

TEST_F(LDTreeTest, BuildBasic) {
  LDTree ldtree(points);

  ldtree.build();

  EXPECT_TRUE(ldtree.tree().isBuilt());
  EXPECT_TRUE(ldtree.voronoi().isBuilt());
}

TEST_F(LDTreeTest, BuildWithVerbose) {
  LDTree ldtree(points);

  std::ostringstream oss;
  ldtree.build(true, // verbose
               &oss  // outputStream
  );

  EXPECT_TRUE(ldtree.tree().isBuilt());
  EXPECT_FALSE(oss.str().empty());
}

TEST_F(LDTreeTest, BuildWithTolerance) {
  LDTree ldtree(points);

  ldtree.build(false,                                   // verbose
               &std::cout,                              // outputStream
               5.0,                                     // logInterval
               std::numeric_limits<double>::infinity(), // timeLimit
               1e-7                                     // tolerance
  );

  EXPECT_TRUE(ldtree.tree().isBuilt());

  EXPECT_THROW(ldtree.build(false, &std::cout, 5.0,
                            std::numeric_limits<double>::infinity(),
                            1e-10 // Too small tolerance
                            ),
               std::invalid_argument);

  EXPECT_THROW(ldtree.build(false, &std::cout, 5.0,
                            std::numeric_limits<double>::infinity(),
                            -1e-1 // Invalid tolerance
                            ),
               std::invalid_argument);
}

TEST_F(LDTreeTest, BuildStats) {
  LDTree ldtree(points);
  ldtree.build();

  const LDTreeStats &stats = ldtree.stats();
  EXPECT_GE(stats.buildTime, 0.0);
}

// ============================================================================
// Query Tests
// ============================================================================

TEST_F(LDTreeTest, QueryBasic) {
  LDTree ldtree(points);
  ldtree.build();

  Index dim = ldtree.voronoi().dimPoints();
  RealVector cost(dim, 1.0);

  std::vector<BinaryVector> result = ldtree.query(cost);

  EXPECT_FALSE(result.empty());

  // Verify result is in the feasible set
  for (const auto &sol : result) {
    bool found = false;
    for (const auto &pt : points) {
      if (sol == pt) {
        found = true;
        break;
      }
    }
    EXPECT_TRUE(found);
  }
}

TEST_F(LDTreeTest, QueryDifferentCosts) {
  LDTree ldtree(points);
  ldtree.build();

  Index dim = ldtree.voronoi().dimPoints();

  RealVector cost1(dim, 1.0);
  RealVector cost2(dim, -1.0);

  std::vector<BinaryVector> result1 = ldtree.query(cost1);
  std::vector<BinaryVector> result2 = ldtree.query(cost2);

  EXPECT_FALSE(result1.empty());
  EXPECT_FALSE(result2.empty());
}

TEST_F(LDTreeTest, QueryResultIsBinary) {
  LDTree ldtree(points);
  ldtree.build();

  Index dim = ldtree.voronoi().dimPoints();
  RealVector cost(dim, 1.0);

  std::vector<BinaryVector> result = ldtree.query(cost);

  for (const auto &sol : result) {
    for (auto val : sol) {
      EXPECT_TRUE(val == 0 || val == 1);
    }
  }
}

// ============================================================================
// Domain Tests
// ============================================================================

TEST_F(LDTreeTest, DomainAccessor) {
  Domain domain;
  domain.push_back(
      {RealVector{1.0, -1.0, 0.0, 0.0, 0.0, 0.0}, 0.0, Relation::GE});

  LDTree ldtree(points, domain);

  EXPECT_EQ(ldtree.domain().size(), 1);
}

// ============================================================================
// Voronoi Accessor Tests
// ============================================================================

TEST_F(LDTreeTest, VoronoiAccessor) {
  LDTree ldtree(points);
  ldtree.build();

  const Voronoi &voronoi = ldtree.voronoi();
  EXPECT_TRUE(voronoi.isBuilt());
  EXPECT_EQ(voronoi.numPoints(), points.size());
}

// ============================================================================
// Tree Accessor Tests
// ============================================================================

TEST_F(LDTreeTest, TreeAccessor) {
  LDTree ldtree(points);
  ldtree.build();

  const Tree &tree = ldtree.tree();
  EXPECT_TRUE(tree.isBuilt());
  EXPECT_GT(tree.size(), 0);
}

// ============================================================================
// Pretty Print Tests
// ============================================================================

TEST_F(LDTreeTest, PprintBasic) {
  LDTree ldtree(points);
  ldtree.build();

  // Should not throw
  EXPECT_NO_THROW(ldtree.pprint(false));
}

TEST_F(LDTreeTest, PprintTight) {
  LDTree ldtree(points);
  ldtree.build();

  // Should not throw
  EXPECT_NO_THROW(ldtree.pprint(true));
}

// ============================================================================
// Flatten Tests
// ============================================================================

TEST_F(LDTreeTest, FlattenBasic) {
  LDTree ldtree(points);
  ldtree.build();

  ldtree.flatten(tempOutputFile);

  // Verify file was created
  std::ifstream infile(tempOutputFile);
  EXPECT_TRUE(infile.good());

  // Verify file contains C code
  std::stringstream buffer;
  buffer << infile.rdbuf();
  std::string content = buffer.str();

  EXPECT_TRUE(content.find("#include") != std::string::npos);
  EXPECT_TRUE(content.find("query") != std::string::npos);
}

TEST_F(LDTreeTest, FlattenWithDoc) {
  LDTree ldtree(points);
  ldtree.build();

  std::string doc = "Test documentation string";
  ldtree.flatten(tempOutputFile, doc);

  std::ifstream infile(tempOutputFile);
  std::stringstream buffer;
  buffer << infile.rdbuf();
  std::string content = buffer.str();

  EXPECT_TRUE(content.find(doc) != std::string::npos);
}

TEST_F(LDTreeTest, FlattenBenchmarkMode) {
  LDTree ldtree(points);
  ldtree.build();

  ldtree.flatten(tempOutputFile, "", true);

  std::ifstream infile(tempOutputFile);
  std::stringstream buffer;
  buffer << infile.rdbuf();
  std::string content = buffer.str();

  // Benchmark mode should include timing code
  EXPECT_TRUE(content.find("time.h") != std::string::npos);
}

// ============================================================================
// Stream Operator Tests
// ============================================================================

TEST_F(LDTreeTest, StreamOperator) {
  LDTree ldtree(points);
  ldtree.build();

  std::ostringstream oss;
  oss << ldtree;
  EXPECT_FALSE(oss.str().empty());
}

// ============================================================================
// Small Instance Tests
// ============================================================================

TEST_F(LDTreeTest, SmallExplicitProblem) {
  std::vector<BinaryVector> smallPoints = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};

  LDTree ldtree(smallPoints);
  ldtree.build();

  EXPECT_TRUE(ldtree.tree().isBuilt());

  RealVector cost = {1.0, 1.0};
  std::vector<BinaryVector> result = ldtree.query(cost);
  EXPECT_FALSE(result.empty());
}

TEST_F(LDTreeTest, TwoPointProblem) {
  std::vector<BinaryVector> twoPoints = {{0, 0, 0}, {1, 1, 1}};

  LDTree ldtree(twoPoints);
  ldtree.build();

  EXPECT_TRUE(ldtree.tree().isBuilt());

  // Query towards {1,1,1}
  RealVector cost1 = {1.0, 1.0, 1.0};
  std::vector<BinaryVector> result1 = ldtree.query(cost1);
  EXPECT_FALSE(result1.empty());
  EXPECT_EQ(result1[0], (BinaryVector{1, 1, 1}));

  // Query towards {0,0,0}
  RealVector cost2 = {-1.0, -1.0, -1.0};
  std::vector<BinaryVector> result2 = ldtree.query(cost2);
  EXPECT_FALSE(result2.empty());
  EXPECT_EQ(result2[0], (BinaryVector{0, 0, 0}));
}
