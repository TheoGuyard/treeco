/**
 * @file Tree.hpp
 * @brief Decision tree structure for optimal solution retrieval.
 *
 * This header defines the Tree class which represents the final decision
 * tree structure that can be queried to find optimal solutions.
 */

#ifndef TREECO_TREE_HPP
#define TREECO_TREE_HPP

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
#include <variant>
#include <vector>

#ifdef TREECO_BUILD_PYTHON
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include "treeco/Dynprog.hpp"
#include "treeco/Geometry.hpp"
#include "treeco/Types.hpp"

namespace treeco {

/// Type of node in the decision tree
enum class NodeType {
  UNDEFINED,  // Not yet determined
  NODE,       // Internal decision node
  LEAF        // Terminal leaf node
};

/**
 * @brief A node in the decision tree.
 *
 * Represents either an internal decision node (with a split and children)
 * or a leaf node (with optimal solution points).
 */
struct Node {
  Index depth = INVALID_INDEX;              // Depth in the tree (0 = root)
  NodeType type = NodeType::UNDEFINED;      // Node type
  std::vector<Index> pointsIds = {};        // Solution point indices (only for LEAF)
  Index splitId = INVALID_INDEX;            // Split hyperplane index (only for NODE)
  std::map<Relation, Index> childIds = {};  // Child node indices by relation (only for NODE)
};

/**
 * @brief Parameters for tree construction.
 */
struct TreeParams {
  bool verbose = false;                                        // Enable verbose logging
  std::ostream* outputStream = &std::cout;                     // Output stream for logs
  double logInterval = 5.0;                                    // Logging interval in seconds
  double timeLimit = std::numeric_limits<double>::infinity();  // Time limit in seconds
  double tolerance = 1e-8;                                     // Numerical tolerance
};

/**
 * @brief Statistics from tree construction.
 */
struct TreeStats {
  bool isBuilt = false;       // Whether the tree has been built
  double buildTime = 0.0;     // Build time in seconds
  DynprogStats dynprogStats;  // Statistics from dynamic programming
  DynprogLogs dynprogLogs;    // Logs from dynamic programming
};

/**
 * @brief Decision tree for querying optimal solutions.
 *
 * This class represents a linear decision tree that can be traversed to
 * find optimal solutions for any given cost vector. The tree is synthesized
 * from the results of dynamic programming over the Voronoi diagram.
 */
class Tree {
public:
  /**
   * @brief Construct a tree for a Voronoi diagram.
   * @param voronoi The Voronoi diagram (must be built)
   */
  Tree(const Voronoi& voronoi);

  /**
   * @brief Synthesize the tree from dynamic programming results.
   * @param dynprog The dynamic programming solver (must have run)
   * @param params Construction parameters
   */
  void synthetize(const Dynprog& dynprog, const TreeParams& params = TreeParams());

  /**
   * @brief Query the tree for optimal solutions.
   * @param cost The cost vector to query
   * @return Vector of optimal solutions (as simplex vectors)
   */
  std::vector<SimplexVector> query(const std::vector<double>& cost) const;

  /**
   * @brief Pretty-print the tree structure.
   * @param tightDisplay If true, show only indices; otherwise show full vectors
   * @param outputStream Output stream for printing
   */
  void pprint(bool tightDisplay = false, std::ostream* outputStream = &std::cout) const;

  /// Clear the tree structure
  void clear();

  /// Check if the tree has been built
  bool isBuilt() const { return stats_.isBuilt; }

  /// Get all nodes
  const std::vector<Node>& nodes() const { return nodes_; }

  /// Get a node by index
  const Node& node(Index id) const { return nodes_.at(id); }

  /// Get the number of nodes
  Index size() const { return size_; }

  /// Get the maximum width (nodes at any level)
  Index width() const { return width_; }

  /// Get the maximum depth
  Index depth() const { return depth_; }

  /// Get construction statistics
  const TreeStats& stats() const { return stats_; }

  /// Get construction parameters
  const TreeParams& params() const { return params_; }

private:
  const Voronoi& voronoi_;   // Reference to the Voronoi diagram
  std::vector<Node> nodes_;  // Tree nodes
  TreeParams params_;        // Construction parameters
  TreeStats stats_;          // Construction statistics

  Index rootId_ = INVALID_INDEX;  // Root node index
  Index size_ = 0;                // Total number of nodes
  Index width_ = 0;               // Maximum width
  Index depth_ = 0;               // Maximum depth

  Clock::time_point startTime_;  // Build start time
  Clock::time_point checkTime_;  // Last logging time

  Index addRoot(const std::vector<Index>& pointsIds, Index splitId);
  Index addNode(Index parentId, Relation relation, const std::vector<Index>& pointsIds, Index splitId);

  void pprintNode(const Node& node, const std::string& prefix, const std::string& label, bool last, bool tight_display,
                  std::ostream* outputStream) const;

  void logHeader() const;
  void logProgress(const std::string& message = "");
  void logFooter() const;
};

/// Stream output operator for Tree
std::ostream& operator<<(std::ostream& oss, const Tree& tree);

};  // namespace treeco

#endif  // TREECO_TREE_HPP