/**
 * @file Voronoi.hpp
 * @brief Voronoi diagram construction for hypercube point sets.
 *
 * This header defines the Voronoi diagram structure used to partition
 * the cost space based on optimality regions of feasible solutions.
 */

#ifndef TREECO_VORONOI_HPP
#define TREECO_VORONOI_HPP

#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#ifdef TREECO_BUILD_PYTHON
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include "treeco/Adjacency.hpp"
#include "treeco/Geometry.hpp"
#include "treeco/Types.hpp"

namespace treeco {

/**
 * @brief A face of the Voronoi diagram.
 *
 * Represents a region of the cost space where a specific feasible solution
 * is optimal. Defined by a point index and the cone defining the face.
 */
struct Face {
  Index pointId = INVALID_INDEX; // Index of the optimal point for this face
  Cone cone = Cone();            // Cone defining the face

  /// Default constructor
  Face() = default;

  /**
   * @brief Construct a face.
   * @param pointId Index of the optimal point
   * @param cone Cone defining the face
   */
  Face(Index pointId, const Cone &cone) : pointId(pointId), cone(cone) {}
};

/**
 * @brief An edge of the Voronoi diagram.
 *
 * Represents the boundary between two adjacent faces, defined by the
 * bisector hyperplane between two optimal points.
 */
struct Edge {
  Index splitId = INVALID_INDEX;  // Index of the split (bisector hyperplane)
  Index leFaceId = INVALID_INDEX; // Face on the LE side of the split
  Index geFaceId = INVALID_INDEX; // Face on the GE side of the split

  /// Default constructor
  Edge() = default;

  /**
   * @brief Construct an edge.
   * @param splitId Index of the split hyperplane
   * @param leFaceId Index of the face on the LE side
   * @param geFaceId Index of the face on the GE side
   */
  Edge(Index splitId, Index leFaceId, Index geFaceId)
      : splitId(splitId), leFaceId(leFaceId), geFaceId(geFaceId) {}
};

/**
 * @brief Parameters for Voronoi diagram construction.
 */
struct VoronoiParams {
  bool verbose = false;                    // Enable verbose logging
  std::ostream *outputStream = &std::cout; // Output stream for logs
  double logInterval = 5.0;                // Logging interval in seconds
  double timeLimit =
      std::numeric_limits<double>::infinity(); // Time limit in seconds
  double tolerance = 1e-8;                     // Numerical tolerance
  bool deduplicate = true;                     // Remove duplicate points
};

/**
 * @brief Statistics from Voronoi diagram construction.
 */
struct VoronoiStats {
  bool isBuilt = false;   // Whether the diagram has been built
  double buildTime = 0.0; // Build time in seconds
  Index lpSolved = 0;     // Number of LP problems solved
};

/**
 * @brief Voronoi diagram for a set of points on the hypercube {-1,+1}^n.
 *
 * This class constructs and maintains the combinatorial structure of the
 * Voronoi diagram, including faces (optimality regions), edges (boundaries),
 * and the relationship between them. The diagram is used to identify which
 * feasible solution is optimal for any given cost vector.
 */
class Voronoi {
public:
  /**
   * @brief Construct a Voronoi diagram for a point set.
   * @param points Set of simplex points in {-1,+1}^n
   */
  Voronoi(const std::vector<SimplexVector> &points);

  /**
   * @brief Build the Voronoi diagram.
   * @param params Construction parameters
   */
  void build(const VoronoiParams &params = VoronoiParams());

  /// Clear the diagram structure and statistics
  void clear();

  /// Check if the diagram has been built
  bool isBuilt() const { return stats_.isBuilt; }

  /// Get the dimension of the points
  Index dimPoints() const { return points_.empty() ? 0 : points_[0].size(); }

  /// Get the number of points
  Index numPoints() const { return points_.size(); }

  /// Get the number of splits (bisector hyperplanes)
  Index numSplits() const { return splits_.size(); }

  /// Get the number of faces
  Index numFaces() const { return faces_.size(); }

  /// Get the number of edges
  Index numEdges() const { return edges_.size(); }

  /// Get a point by index
  const SimplexVector &point(Index id) const { return points_[id]; }

  /// Get a split by index
  const TernaryVector &split(Index id) const { return splits_[id]; }

  /// Get a face by index
  const Face &face(Index id) const { return faces_[id]; }

  /// Get an edge by index
  const Edge &edge(Index id) const { return edges_[id]; }

  /// Get all points
  const std::vector<SimplexVector> &points() const { return points_; }

  /// Get all splits
  const std::vector<TernaryVector> &splits() const { return splits_; }

  /// Get all faces
  const std::vector<Face> &faces() const { return faces_; }

  /// Get all edges
  const std::vector<Edge> &edges() const { return edges_; }

  /// Get construction statistics
  const VoronoiStats &stats() const { return stats_; }

  /// Get construction parameters
  const VoronoiParams &params() const { return params_; }

private:
  const std::vector<SimplexVector> points_; // Input points
  std::vector<TernaryVector> splits_;       // Bisector hyperplanes
  std::vector<Face> faces_;                 // Voronoi faces
  std::vector<Edge> edges_;                 // Voronoi edges
  VoronoiParams params_;                    // Construction parameters
  VoronoiStats stats_;                      // Construction statistics

  Clock::time_point startTime_; // Build start time
  Clock::time_point checkTime_; // Last logging time

  void logHeader() const;
  void logProgress(const Adjacency &adj, Index i,
                   const std::string &message = "");
  void logFooter() const;
};

/// Stream output operator for Voronoi diagram
std::ostream &operator<<(std::ostream &oss, const Voronoi &voronoi);

} // namespace treeco

#endif // TREECO_VORONOI_HPP
