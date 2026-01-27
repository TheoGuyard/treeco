#include <gtest/gtest.h>
#include <random>
#include <treeco/Geometry.hpp>
#include <treeco/Problem/Maxcut.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Types.hpp>
#include <treeco/Voronoi.hpp>

using namespace treeco;

class VoronoiTest : public ::testing::Test {
protected:
  std::vector<SimplexVector> points;

  void SetUp() override {
    Tsp tsp(5);
    points = scaleBinarySet(tsp.getFeasibleSet());
  }
};

TEST_F(VoronoiTest, BasicConstruction) {
  Voronoi voronoi(points);
  voronoi.build();

  EXPECT_EQ(voronoi.numFaces(), voronoi.numPoints());

  Index totalCuts = 0;
  for (Index i = 0; i < voronoi.numFaces(); ++i) {
    totalCuts += voronoi.face(i).cone.numCuts();
  }
  EXPECT_EQ(voronoi.numEdges(), totalCuts / 2);

  for (Index i = 0; i < voronoi.numSplits(); ++i) {
    auto split = voronoi.split(i);
    for (auto val : split) {
      if (val != 0) {
        EXPECT_EQ(val, 1);
        break;
      }
    }
  }

  VoronoiParams params;
  params.deduplicate = false;
  Voronoi voronoiDuplicate(points);
  voronoiDuplicate.build(params);

  EXPECT_GE(voronoiDuplicate.numSplits(), voronoi.numSplits());
}

TEST_F(VoronoiTest, CellsMatchPoints) {
  Voronoi voronoi(points);
  voronoi.build();

  for (Index i = 0; i < voronoi.numFaces(); ++i) {
    EXPECT_EQ(voronoi.face(i).pointId, i);
  }
}
