#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <treeco/IO.hpp>
#include <treeco/Types.hpp>

using namespace treeco;

class IOTest : public ::testing::Test {
protected:
  std::string tempFilePath;

  void SetUp() override { tempFilePath = "/tmp/treeco_test_io_" + std::to_string(std::rand()) + ".txt"; }

  void TearDown() override { std::remove(tempFilePath.c_str()); }

  void writeToFile(const std::string& content) {
    std::ofstream outfile(tempFilePath);
    outfile << content;
    outfile.close();
  }
};

// ============================================================================
// readPoints Tests
// ============================================================================

TEST_F(IOTest, ReadPointsBasic) {
  writeToFile("points 3 4\n"
              "0 0 0\n"
              "1 0 0\n"
              "0 1 0\n"
              "1 1 1\n");

  std::vector<BinaryVector> points = readPoints(tempFilePath);

  ASSERT_EQ(points.size(), 4);
  ASSERT_EQ(points[0].size(), 3);

  EXPECT_EQ(points[0], (BinaryVector{0, 0, 0}));
  EXPECT_EQ(points[1], (BinaryVector{1, 0, 0}));
  EXPECT_EQ(points[2], (BinaryVector{0, 1, 0}));
  EXPECT_EQ(points[3], (BinaryVector{1, 1, 1}));
}

TEST_F(IOTest, ReadPointsSinglePoint) {
  writeToFile("points 2 1\n"
              "1 0\n");

  std::vector<BinaryVector> points = readPoints(tempFilePath);

  ASSERT_EQ(points.size(), 1);
  EXPECT_EQ(points[0], (BinaryVector{1, 0}));
}

TEST_F(IOTest, ReadPointsWithEmptyLines) {
  writeToFile("points 2 2\n"
              "\n"
              "0 1\n"
              "\n"
              "1 0\n"
              "\n");

  std::vector<BinaryVector> points = readPoints(tempFilePath);

  ASSERT_EQ(points.size(), 2);
  EXPECT_EQ(points[0], (BinaryVector{0, 1}));
  EXPECT_EQ(points[1], (BinaryVector{1, 0}));
}

TEST_F(IOTest, ReadPointsFileNotFound) { EXPECT_THROW(readPoints("/nonexistent/path/file.txt"), std::runtime_error); }

TEST_F(IOTest, ReadPointsEmptyFile) {
  writeToFile("");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsInvalidHeader) {
  writeToFile("invalid header format\n0 0\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsWrongKeyword) {
  writeToFile("vertices 2 1\n0 1\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsZeroDimension) {
  writeToFile("points 0 1\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsZeroCount) {
  writeToFile("points 2 0\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsNonBinaryEntry) {
  writeToFile("points 2 1\n"
              "0 2\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsWrongDimension) {
  writeToFile("points 3 1\n"
              "0 1\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadPointsTooFewPoints) {
  writeToFile("points 2 3\n"
              "0 1\n"
              "1 0\n");
  EXPECT_THROW(readPoints(tempFilePath), std::runtime_error);
}

// ============================================================================
// writePoints Tests
// ============================================================================

TEST_F(IOTest, WritePointsBasic) {
  std::vector<BinaryVector> points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 1}};

  writePoints(tempFilePath, points);

  // Read back and verify
  std::vector<BinaryVector> readBack = readPoints(tempFilePath);

  ASSERT_EQ(readBack.size(), 4);
  EXPECT_EQ(readBack[0], (BinaryVector{0, 0, 0}));
  EXPECT_EQ(readBack[1], (BinaryVector{1, 0, 0}));
  EXPECT_EQ(readBack[2], (BinaryVector{0, 1, 0}));
  EXPECT_EQ(readBack[3], (BinaryVector{1, 1, 1}));
}

TEST_F(IOTest, WritePointsEmpty) {
  std::vector<BinaryVector> points;
  EXPECT_THROW(writePoints(tempFilePath, points), std::runtime_error);
}

TEST_F(IOTest, WritePointsInconsistentDimensions) {
  std::vector<BinaryVector> points = {
      {0, 0, 0}, {1, 0}  // Different dimension
  };
  EXPECT_THROW(writePoints(tempFilePath, points), std::runtime_error);
}

TEST_F(IOTest, WriteReadPointsRoundTrip) {
  std::vector<BinaryVector> original = {{1, 0, 1, 0, 1}, {0, 1, 0, 1, 0}, {1, 1, 1, 1, 1}};

  writePoints(tempFilePath, original);
  std::vector<BinaryVector> readBack = readPoints(tempFilePath);

  ASSERT_EQ(readBack.size(), original.size());
  for (size_t i = 0; i < original.size(); ++i) { EXPECT_EQ(readBack[i], original[i]); }
}

// ============================================================================
// readDomain Tests
// ============================================================================

TEST_F(IOTest, ReadDomainBasic) {
  writeToFile("domain 3 2\n"
              "1.0 2.0 -1.0 3.5 LT\n"
              "0.5 -0.5 1.0 0.0 EQ\n");

  Domain domain = readDomain(tempFilePath);

  ASSERT_EQ(domain.size(), 2);

  // First constraint
  const auto& [coeffs1, b1, rel1] = domain[0];
  ASSERT_EQ(coeffs1.size(), 3);
  EXPECT_DOUBLE_EQ(coeffs1[0], 1.0);
  EXPECT_DOUBLE_EQ(coeffs1[1], 2.0);
  EXPECT_DOUBLE_EQ(coeffs1[2], -1.0);
  EXPECT_DOUBLE_EQ(b1, 3.5);
  EXPECT_EQ(rel1, Relation::LT);

  // Second constraint
  const auto& [coeffs2, b2, rel2] = domain[1];
  ASSERT_EQ(coeffs2.size(), 3);
  EXPECT_DOUBLE_EQ(coeffs2[0], 0.5);
  EXPECT_DOUBLE_EQ(coeffs2[1], -0.5);
  EXPECT_DOUBLE_EQ(coeffs2[2], 1.0);
  EXPECT_DOUBLE_EQ(b2, 0.0);
  EXPECT_EQ(rel2, Relation::EQ);
}

TEST_F(IOTest, ReadDomainAllRelations) {
  writeToFile("domain 2 5\n"
              "1.0 0.0 0.0 LT\n"
              "0.0 1.0 0.0 LE\n"
              "1.0 1.0 0.0 EQ\n"
              "-1.0 0.0 0.0 GE\n"
              "0.0 -1.0 0.0 GT\n");

  Domain domain = readDomain(tempFilePath);

  ASSERT_EQ(domain.size(), 5);
  EXPECT_EQ(std::get<2>(domain[0]), Relation::LT);
  EXPECT_EQ(std::get<2>(domain[1]), Relation::LE);
  EXPECT_EQ(std::get<2>(domain[2]), Relation::EQ);
  EXPECT_EQ(std::get<2>(domain[3]), Relation::GE);
  EXPECT_EQ(std::get<2>(domain[4]), Relation::GT);
}

TEST_F(IOTest, ReadDomainFileNotFound) { EXPECT_THROW(readDomain("/nonexistent/path/file.txt"), std::runtime_error); }

TEST_F(IOTest, ReadDomainEmptyFile) {
  writeToFile("");
  EXPECT_THROW(readDomain(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadDomainWrongKeyword) {
  writeToFile("constraints 2 1\n1.0 0.0 0.0 LT\n");
  EXPECT_THROW(readDomain(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadDomainInvalidRelation) {
  writeToFile("domain 2 1\n"
              "1.0 0.0 0.0 ?\n");
  EXPECT_THROW(readDomain(tempFilePath), std::invalid_argument);
}

TEST_F(IOTest, ReadDomainMissingCoefficients) {
  writeToFile("domain 3 1\n"
              "1.0 2.0 0.0 LT\n"  // Missing b
  );
  EXPECT_THROW(readDomain(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadDomainWrongConstraintCount) {
  writeToFile("domain 2 3\n"
              "1.0 0.0 0.0 LT\n"
              "0.0 1.0 0.0 GT\n");
  EXPECT_THROW(readDomain(tempFilePath), std::runtime_error);
}

// ============================================================================
// writeDomain Tests
// ============================================================================

TEST_F(IOTest, WriteDomainBasic) {
  Domain domain = {{RealVector{1.0, 2.0, -1.0}, 3.5, Relation::LT}, {RealVector{0.5, -0.5, 1.0}, 0.0, Relation::EQ}};

  writeDomain(tempFilePath, domain);

  // Read back and verify
  Domain readBack = readDomain(tempFilePath);

  ASSERT_EQ(readBack.size(), 2);

  const auto& [coeffs1, b1, rel1] = readBack[0];
  EXPECT_DOUBLE_EQ(coeffs1[0], 1.0);
  EXPECT_DOUBLE_EQ(coeffs1[1], 2.0);
  EXPECT_DOUBLE_EQ(coeffs1[2], -1.0);
  EXPECT_DOUBLE_EQ(b1, 3.5);
  EXPECT_EQ(rel1, Relation::LT);
}

TEST_F(IOTest, WriteDomainAllRelations) {
  Domain domain = {{RealVector{1.0, 0.0}, 0.0, Relation::LT},
                   {RealVector{0.0, 1.0}, 0.0, Relation::LE},
                   {RealVector{1.0, 1.0}, 0.0, Relation::EQ},
                   {RealVector{-1.0, 0.0}, 0.0, Relation::GE},
                   {RealVector{0.0, -1.0}, 0.0, Relation::GT}};

  writeDomain(tempFilePath, domain);
  Domain readBack = readDomain(tempFilePath);

  ASSERT_EQ(readBack.size(), 5);
  EXPECT_EQ(std::get<2>(readBack[0]), Relation::LT);
  EXPECT_EQ(std::get<2>(readBack[1]), Relation::LE);
  EXPECT_EQ(std::get<2>(readBack[2]), Relation::EQ);
  EXPECT_EQ(std::get<2>(readBack[3]), Relation::GE);
  EXPECT_EQ(std::get<2>(readBack[4]), Relation::GT);
}

TEST_F(IOTest, WriteDomainEmpty) {
  Domain domain;
  EXPECT_THROW(writeDomain(tempFilePath, domain), std::runtime_error);
}

TEST_F(IOTest, WriteDomainInconsistentDimensions) {
  Domain domain = {
      {RealVector{1.0, 2.0, 3.0}, 0.0, Relation::LT}, {RealVector{1.0, 2.0}, 0.0, Relation::GT}  // Different dimension
  };
  EXPECT_THROW(writeDomain(tempFilePath, domain), std::runtime_error);
}

// ============================================================================
// readQueries Tests
// ============================================================================

TEST_F(IOTest, ReadQueriesBasic) {
  writeToFile("queries 3 2\n"
              "1.0 2.5 -3.0\n"
              "-1.5 0.0 4.5\n");

  std::vector<RealVector> queries = readQueries(tempFilePath);

  ASSERT_EQ(queries.size(), 2);
  ASSERT_EQ(queries[0].size(), 3);

  EXPECT_DOUBLE_EQ(queries[0][0], 1.0);
  EXPECT_DOUBLE_EQ(queries[0][1], 2.5);
  EXPECT_DOUBLE_EQ(queries[0][2], -3.0);
  EXPECT_DOUBLE_EQ(queries[1][0], -1.5);
  EXPECT_DOUBLE_EQ(queries[1][1], 0.0);
  EXPECT_DOUBLE_EQ(queries[1][2], 4.5);
}

TEST_F(IOTest, ReadQueriesSingleQuery) {
  writeToFile("queries 2 1\n"
              "1.5 -2.5\n");

  std::vector<RealVector> queries = readQueries(tempFilePath);

  ASSERT_EQ(queries.size(), 1);
  EXPECT_DOUBLE_EQ(queries[0][0], 1.5);
  EXPECT_DOUBLE_EQ(queries[0][1], -2.5);
}

TEST_F(IOTest, ReadQueriesWithEmptyLines) {
  writeToFile("queries 2 2\n"
              "\n"
              "1.0 2.0\n"
              "\n"
              "3.0 4.0\n"
              "\n");

  std::vector<RealVector> queries = readQueries(tempFilePath);

  ASSERT_EQ(queries.size(), 2);
  EXPECT_DOUBLE_EQ(queries[0][0], 1.0);
  EXPECT_DOUBLE_EQ(queries[1][0], 3.0);
}

TEST_F(IOTest, ReadQueriesFileNotFound) { EXPECT_THROW(readQueries("/nonexistent/path/file.txt"), std::runtime_error); }

TEST_F(IOTest, ReadQueriesEmptyFile) {
  writeToFile("");
  EXPECT_THROW(readQueries(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadQueriesWrongKeyword) {
  writeToFile("costs 2 1\n1.0 2.0\n");
  EXPECT_THROW(readQueries(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadQueriesZeroDimension) {
  writeToFile("queries 0 1\n");
  EXPECT_THROW(readQueries(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadQueriesWrongDimension) {
  writeToFile("queries 3 1\n"
              "1.0 2.0\n");
  EXPECT_THROW(readQueries(tempFilePath), std::runtime_error);
}

TEST_F(IOTest, ReadQueriesTooFewQueries) {
  writeToFile("queries 2 3\n"
              "1.0 2.0\n"
              "3.0 4.0\n");
  EXPECT_THROW(readQueries(tempFilePath), std::runtime_error);
}

// ============================================================================
// writeQueries Tests
// ============================================================================

TEST_F(IOTest, WriteQueriesBasic) {
  std::vector<RealVector> queries = {{1.0, 2.5, -3.0}, {-1.5, 0.0, 4.5}};

  writeQueries(tempFilePath, queries);

  // Read back and verify
  std::vector<RealVector> readBack = readQueries(tempFilePath);

  ASSERT_EQ(readBack.size(), 2);
  EXPECT_DOUBLE_EQ(readBack[0][0], 1.0);
  EXPECT_DOUBLE_EQ(readBack[0][1], 2.5);
  EXPECT_DOUBLE_EQ(readBack[0][2], -3.0);
  EXPECT_DOUBLE_EQ(readBack[1][0], -1.5);
  EXPECT_DOUBLE_EQ(readBack[1][1], 0.0);
  EXPECT_DOUBLE_EQ(readBack[1][2], 4.5);
}

TEST_F(IOTest, WriteQueriesEmpty) {
  std::vector<RealVector> queries;
  EXPECT_THROW(writeQueries(tempFilePath, queries), std::runtime_error);
}

TEST_F(IOTest, WriteQueriesInconsistentDimensions) {
  std::vector<RealVector> queries = {
      {1.0, 2.0, 3.0}, {1.0, 2.0}  // Different dimension
  };
  EXPECT_THROW(writeQueries(tempFilePath, queries), std::runtime_error);
}

TEST_F(IOTest, WriteReadQueriesRoundTrip) {
  std::vector<RealVector> original = {{1.23456789, -9.87654321, 0.0}, {3.14159265, 2.71828182, -1.41421356}};

  writeQueries(tempFilePath, original);
  std::vector<RealVector> readBack = readQueries(tempFilePath);

  ASSERT_EQ(readBack.size(), original.size());
  for (size_t i = 0; i < original.size(); ++i) {
    ASSERT_EQ(readBack[i].size(), original[i].size());
    for (size_t j = 0; j < original[i].size(); ++j) { EXPECT_NEAR(readBack[i][j], original[i][j], 1e-10); }
  }
}

TEST_F(IOTest, WriteReadQueriesPrecision) {
  // Test that high precision values are preserved
  std::vector<RealVector> original = {{1.234567890123456, -9.876543210987654}};

  writeQueries(tempFilePath, original);
  std::vector<RealVector> readBack = readQueries(tempFilePath);

  ASSERT_EQ(readBack.size(), 1);
  EXPECT_NEAR(readBack[0][0], original[0][0], 1e-15);
  EXPECT_NEAR(readBack[0][1], original[0][1], 1e-15);
}