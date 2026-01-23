#include <gtest/gtest.h>
#include <treeco/Adjacency.hpp>
#include <treeco/Geometry.hpp>
#include <treeco/Problem/Tsp.hpp>

using namespace treeco;

class AdjacencyTest : public ::testing::Test {
protected:

    std::vector<SimplexVector> points;
    
    void SetUp() override {
        Tsp tsp(5);
        points = scaleBinarySet(tsp.getFeasibleSet());
    }
};


TEST_F(AdjacencyTest, AdjacencyChecks) {
    Adjacency adj(points);  
    EXPECT_TRUE(adj.check(0, 1));
    EXPECT_TRUE(adj.check(0, 2));
    EXPECT_TRUE(adj.check(1, 2));
}

TEST_F(AdjacencyTest, SelfAdjacentError) {
    Adjacency adj(points);  
    EXPECT_THROW(adj.check(0, 0), std::invalid_argument);
    EXPECT_THROW(adj.check(1, 1), std::invalid_argument);
}
