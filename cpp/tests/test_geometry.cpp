#include <gtest/gtest.h>

#include <treeco/Geometry.hpp>

using namespace treeco;

TEST(GeometryTest, ScaleBinary) {
  BinaryVector x = {0, 1, 0, 1};
  SimplexVector y = scaleBinary(x);

  ASSERT_EQ(y.size(), 4);
  EXPECT_EQ(y[0], -1);
  EXPECT_EQ(y[1], 1);
  EXPECT_EQ(y[2], -1);
  EXPECT_EQ(y[3], 1);
}

TEST(GeometryTest, TernaryNormalization) {
  TernaryVector v1 = {0, -1, 1};
  normalize(v1);
  EXPECT_EQ(v1[0], 0);
  EXPECT_EQ(v1[1], 1);
  EXPECT_EQ(v1[2], -1);

  TernaryVector v2 = {0, 1, -1};
  normalize(v2);
  EXPECT_EQ(v2[0], 0);
  EXPECT_EQ(v2[1], 1);
  EXPECT_EQ(v2[2], -1);

  TernaryVector v3 = {1, 0, 1};
  normalize(v3);
  EXPECT_EQ(v3[0], 1);
  EXPECT_EQ(v3[1], 0);
  EXPECT_EQ(v3[2], 1);
}

TEST(GeometryTest, Bisector) {
  SimplexVector p1 = {-1, -1, 1};
  SimplexVector p2 = {1, -1, 1};

  TernaryVector v = bisector(p1, p2);

  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 0);
  EXPECT_EQ(v[2], 0);
}

TEST(GeometryTest, ConeHash) {
  Cone c;
  c.addCut(0, Relation::LT);
  c.addCut(1, Relation::GT);

  std::size_t hash = c.hash();
  EXPECT_NE(hash, 0);

  Cone c2;
  c2.addCut(0, Relation::LT);
  c2.addCut(1, Relation::GT);

  EXPECT_EQ(c.hash(), c2.hash());

  Cone c3;
  c3.addCut(0, Relation::GT);
  c3.addCut(1, Relation::LT);

  EXPECT_NE(c.hash(), c3.hash());
}
