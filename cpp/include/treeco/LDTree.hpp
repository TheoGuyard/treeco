/**
 * @file LDTree.hpp
 * @brief Linear Decision Tree for combinatorial optimization.
 *
 * This is the main user-facing class that provides a high-level interface
 * for building and querying linear decision trees for optimization problems.
 */

#ifndef TREECO_LDTREE_HPP
#define TREECO_LDTREE_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "treeco/Dynprog.hpp"
#include "treeco/Geometry.hpp"
#include "treeco/IO.hpp"
#include "treeco/Tree.hpp"
#include "treeco/Types.hpp"
#include "treeco/Voronoi.hpp"

namespace treeco {

/**
 * @brief Statistics from LDTree construction.
 */
struct LDTreeStats {
  double buildTime = 0.0;  // Total build time in seconds
};

/**
 * @brief Linear Decision Tree for combinatorial optimization.
 *
 * This class provides a complete workflow for building decision trees that
 * can efficiently solve linear optimization problems over a fixed feasible set.
 * Once built, the tree can be queried with any cost vector to retrieve the
 * optimal solutions in O(depth) time.
 *
 * Example usage:
 * @code
 * LDTree tree(points, domain);
 * tree.build();
 * auto solutions = tree.query(cost);
 * @endcode
 */
class LDTree {
public:
  /**
   * @brief Construct an LDTree from input files.
   * @param filePoints Path to file containing feasible points
   * @param fileDomain Path to file containing cost domain constraints
   * (optional)
   */
  LDTree(const std::string& filePoints, const std::string& fileDomain = "");

  /**
   * @brief Construct an LDTree from data objects.
   * @param points Vector of feasible binary solutions
   * @param domain Optional linear constraints on the cost domain
   */
  LDTree(const std::vector<BinaryVector>& points, const Domain& domain = Domain());

  /// Serialization constructor
  LDTree(const Domain& domain, const LDTreeStats& stats, const std::vector<SimplexVector>& voronoiPoints,
         const std::vector<TernaryVector>& voronoiSplits, const std::vector<Face>& voronoiFaces,
         const std::vector<Edge>& voronoiEdges, const VoronoiParams& voronoiParams, const VoronoiStats& voronoiStats,
         const std::vector<Node>& treeNodes, const TreeParams& treeParams, const TreeStats& treeStats, Index treeRootId,
         Index treeSize, Index treeWidth, Index treeDepth)
    : domain_(domain)
    , voronoi_(voronoiPoints, voronoiSplits, voronoiFaces, voronoiEdges, voronoiParams, voronoiStats)
    , tree_(treeNodes, treeParams, treeStats, treeRootId, treeSize, treeWidth, treeDepth)
    , stats_(stats) {}

  /**
   * @brief Build the LDTree structure.
   *
   * Constructs the Voronoi diagram and decision tree structure using
   * dynamic programming to find the minimum-depth tree.
   *
   * @param verbose Enable verbose logging output
   * @param outputStream Output stream for logging (default: std::cout)
   * @param logInterval Logging interval in seconds
   * @param logSave Save logs at each iteration of the dynamic programming
   * @param timeLimit Maximum build time in seconds
   * @param tolerance Numerical tolerance for equality comparisons
   * @param deduplicate Remove duplicate feasible points before building
   * @param filterChecks Use filtering for faster split validity checks
   * @param exploration Search exploration strategy
   * @param branching Tree branching mode (binary or ternary)
   * @param lowerBounding Lower bound computation strategy
   * @param positioning Face-split position computation mode
   * @param splitSelection Split candidate selection strategy
   * @param splitScoring Split quality scoring strategy
   * @param randomSeed Random seed for sampling-based selection
   */
  void build(bool verbose = false, std::ostream* outputStream = &std::cout, double logInterval = 5.0,
             bool logSave = true, double timeLimit = std::numeric_limits<double>::infinity(), double tolerance = 1e-8,
             bool deduplicate = true, bool filterChecks = true, Exploration exploration = Exploration::ITERATIVE,
             Branching branching = Branching::BINARY, LowerBounding lowerBounding = LowerBounding::BACKTRACK,
             Positioning positioning = Positioning::ONLINE, SplitSelection splitSelection = SplitSelection::ALL,
             SplitScoring splitScoring = SplitScoring::VARIANCE, Index randomSeed = 42);

  /**
   * @brief Query the tree for optimal solutions.
   * @param cost The cost vector to optimize
   * @param tolerance Numerical tolerance for split evaluations
   * @param checkDomain Validate that cost is within the domain before querying
   * @return Vector of optimal binary solutions
   */
  std::vector<BinaryVector> query(const RealVector& cost, double tolerance = 1e-8, bool checkDomain = false) const;

  /**
   * @brief Pretty-print the tree structure.
   * @param tightDisplay If true, show only indices; otherwise show full vectors
   * @param outputStream Output stream for printing
   */
  void pprint(bool tightDisplay = false, std::ostream* outputStream = &std::cout) const;

  /**
   * @brief Generate standalone C code implementing the decision tree.
   * @param filepath Output file path for the generated code
   * @param doc Documentation header for the generated code
   * @param benchmarkMode Generate benchmarking code instead of query code
   */
  void flatten(const std::string filepath, const std::string doc = "", bool benchmarkMode = false) const;

  /// Get the cost domain constraints
  const Domain& domain() const noexcept { return domain_; }

  /// Get the Voronoi diagram
  const Voronoi& voronoi() const noexcept { return voronoi_; }

  /// Get the decision tree structure
  const Tree& tree() const noexcept { return tree_; }

  /// Get build statistics
  const LDTreeStats& stats() const noexcept { return stats_; }

private:
  Domain domain_;      // Cost domain constraints
  Voronoi voronoi_;    // Voronoi diagram of feasible points
  Tree tree_;          // Decision tree structure
  LDTreeStats stats_;  // Build statistics

  // LDTree::pprint helpers
  void pprintNode(const Node& node, const std::string& prefix, const std::string& label, bool last, bool tightDisplay,
                  std::ostream* outputStream) const;

  // LDTree::flatten helpers
  void generateNormalCode(std::ostringstream& out, const std::string& doc) const;
  void generateBenchmarkCode(std::ostringstream& out, const std::string& doc) const;
  void generateNodeCode(Index nodeIndex, std::ostringstream& out, int indent = 1) const;
  void generateDotCode(std::ostringstream& out, const TernaryVector& split) const;
};

/// Stream output operator for LDTree
std::ostream& operator<<(std::ostream& oss, const LDTree& tree);

}  // namespace treeco

#endif  // TREECO_LDTREE_HPP