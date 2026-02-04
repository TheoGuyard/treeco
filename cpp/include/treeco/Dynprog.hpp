/**
 * @file Dynprog.hpp
 * @brief Dynamic programming solver for minimum-depth decision trees.
 *
 * This header defines the dynamic programming algorithm that computes
 * the optimal (minimum-depth) decision tree structure over the Voronoi diagram.
 */

#ifndef TREECO_DYNPROG_HPP
#define TREECO_DYNPROG_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef TREECO_BUILD_PYTHON
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include "treeco/Feasibility.hpp"
#include "treeco/Geometry.hpp"
#include "treeco/Types.hpp"
#include "treeco/Voronoi.hpp"

namespace treeco {

/// Maximum possible depth value (used as infinity)
constexpr Index MAX_DEPTH = std::numeric_limits<Index>::max();

/**
 * @brief Safely add to a depth value, avoiding overflow.
 * @param d Current depth
 * @param increment Value to add (default: 1)
 * @return Sum or MAX_DEPTH if overflow would occur
 */
inline Index safeAdd(Index d, Index increment = 1) {
  return (d >= MAX_DEPTH - increment) ? MAX_DEPTH : d + increment;
}

/**
 * @brief Exploration strategy for the dynamic programming algorithm.
 */
enum class Exploration {
  GREEDY,     // Single iteration with k=1 (fast, potentially suboptimal)
  ITERATIVE,  // Iterative scheme with k=1,2,...,numSplits (optimal, early
              // stopping)
  EXHAUSTIVE  // Single iteration with k=numSplits (optimal, no early stopping)
};

/**
 * @brief Convert Exploration enum to string.
 * @param type The exploration type
 * @return String representation
 */
inline std::string explorationTypeToString(Exploration type) {
  switch (type) {
    case Exploration::GREEDY:
      return "greedy";
    case Exploration::ITERATIVE:
      return "iterative";
    case Exploration::EXHAUSTIVE:
      return "exhaustive";
    default:
      return "unknown";
  }
}

/**
 * @brief Parse string to Exploration enum.
 * @param str String representation
 * @return Corresponding Exploration value
 * @throws std::invalid_argument if string is invalid
 */
inline Exploration stringToExplorationType(const std::string& str) {
  if (str == "greedy") return Exploration::GREEDY;
  if (str == "iterative") return Exploration::ITERATIVE;
  if (str == "exhaustive") return Exploration::EXHAUSTIVE;
  throw std::invalid_argument("Invalid exploration type: " + str +
                              ". Available options are: "
                              "greedy, iterative, exhaustive.");
}

/**
 * @brief Branching mode for decision tree nodes.
 */
enum class Branching {
  TERNARY,  // Three-way branching: {<, =, >}
  BINARY    // Two-way branching: {<, >} with tie-breaking for {=}
};

/**
 * @brief Convert Branching enum to string.
 * @param type The branching type
 * @return String representation
 */
inline std::string branchingTypeToString(Branching type) {
  switch (type) {
    case Branching::BINARY:
      return "binary";
    case Branching::TERNARY:
      return "ternary";
    default:
      return "unknown";
  }
}

/**
 * @brief Parse string to Branching enum.
 * @param str String representation
 * @return Corresponding Branching value
 * @throws std::invalid_argument if string is invalid
 */
inline Branching stringToBranchingType(const std::string& str) {
  if (str == "binary") return Branching::BINARY;
  if (str == "ternary") return Branching::TERNARY;
  throw std::invalid_argument("Invalid branching type: " + str +
                              ". Available options are: "
                              "binary, ternary.");
}

/**
 * @brief Split scoring strategy for ordering candidate splits.
 */
enum class SplitScoring {
  VARIANCE,  // Variance-based deviation from equal split
  ENTROPY,   // Information gain (entropy reduction)
  MINMAX,    // Minimize maximum number of child faces
  NONE,      // All splits have equal score
  RANDOM     // Random scoring
};

/**
 * @brief Convert SplitScoring enum to string.
 * @param type The split scoring type
 * @return String representation
 */
inline std::string splitScoringTypeToString(SplitScoring type) {
  switch (type) {
    case SplitScoring::VARIANCE:
      return "variance";
    case SplitScoring::ENTROPY:
      return "entropy";
    case SplitScoring::MINMAX:
      return "minmax";
    case SplitScoring::NONE:
      return "none";
    case SplitScoring::RANDOM:
      return "random";
    default:
      return "unknown";
  }
}

/**
 * @brief Parse string to SplitScoring enum.
 * @param str String representation
 * @return Corresponding SplitScoring value
 * @throws std::invalid_argument if string is invalid
 */
inline SplitScoring stringToSplitScoringType(const std::string& str) {
  if (str == "variance") return SplitScoring::VARIANCE;
  if (str == "entropy") return SplitScoring::ENTROPY;
  if (str == "minmax") return SplitScoring::MINMAX;
  if (str == "none") return SplitScoring::NONE;
  if (str == "random") return SplitScoring::RANDOM;
  throw std::invalid_argument("Invalid split scoring type: " + str +
                              ". Available options are: "
                              "variance, entropy, minmax, none, random.");
}

/**
 * @brief Statistics from the dynamic programming algorithm.
 */
struct DynprogStats {
  double runTime = 0.0;            // Total running time in seconds
  Index numIters = 0;              // Iterations in iterative scheme
  Index numEvals = 0;              // State evaluations performed
  Index numStates = 0;             // Total states created
  Index numStatesBuilt = 0;        // States fully built
  Index numStatesClosed = 0;       // States closed by optimality
  Index numStatesLeafed = 0;       // States identified as leaves
  Index numStatesPruned = 0;       // States pruned by bounds
  Index lpSolved = 0;              // LP feasibility checks performed
  Index optimalDepth = MAX_DEPTH;  // Optimal tree depth found
};

/// Stream output operator for DynprogStats
std::ostream& operator<<(std::ostream& oss, const DynprogStats& stats);

/// Log entry for dynamic programming progress (sequence of stats)
using DynprogLogs = std::vector<DynprogStats>;

/**
 * @brief Termination status of the dynamic programming algorithm.
 */
enum class DynprogStatus {
  INVALID,     // Not started, in progress, or terminated without valid tree
  SUBOPTIMAL,  // Found valid tree but optimality not certified
  OPTIMAL      // Found and certified optimal tree depth
};

/**
 * @brief Convert DynprogStatus enum to string.
 * @param status The status
 * @return String representation
 */
inline std::string dynprogStatusToString(DynprogStatus status) {
  switch (status) {
    case DynprogStatus::INVALID:
      return "invalid";
    case DynprogStatus::SUBOPTIMAL:
      return "suboptimal";
    case DynprogStatus::OPTIMAL:
      return "optimal";
    default:
      return "unknown";
  }
}

// ----------------------------------------------------------------------------
// Dynamic programming state
// ----------------------------------------------------------------------------

/**
 * @brief Information about a candidate split at a state.
 */
struct Split {
  Index id = INVALID_INDEX;  // Split index in Voronoi diagram
  std::optional<bool> valid =
      std::nullopt;  // Whether the split is valid for this state
  double score =
      std::numeric_limits<double>::infinity();  // Split quality score
  std::array<Index, 3> childIds = {};           // Child state indices
  bool isClosed = false;  // Whether the split is fully explored
};

/**
 * @brief A state in the dynamic programming search.
 *
 * Represents a subproblem defined by a polyhedral region and the set of
 * Voronoi faces that intersect it.
 */
struct State {
  Cone region;                      // Polyhedral region of the state
  std::vector<Index> faceIds = {};  // Voronoi face indices in this region
  std::vector<Split> splits = {};   // Valid splits for this state
  Index lbHeight = 0;               // Lower bound on optimal subtree height
  Index ubHeight = MAX_DEPTH;       // Upper bound on optimal subtree height
  Index splitId = INVALID_INDEX;    // Best split index found
  bool isBuilt = false;             // Whether all splits have been explored
  bool isClosed = false;  // Whether state is closed (leaf, pruned, or optimal)

  /// Get the depth of this state in the tree
  Index depth() const { return static_cast<Index>(region.numCuts()); }

  /// Check if this state is a leaf (single face or no valid splits)
  bool isLeaf() const { return faceIds.size() <= 1 || splits.size() == 0; }
};

// ----------------------------------------------------------------------------
// Dynamic programming method
// ----------------------------------------------------------------------------

/**
 * @brief Parameters for the dynamic programming algorithm.
 */
struct DynprogParams {
  bool verbose = false;                     // Enable verbose logging
  std::ostream* outputStream = &std::cout;  // Output stream for logs
  double logInterval = 5.0;                 // Logging interval in seconds
  bool logSave = false;                     // Save logs at each iteration
  double timeLimit =
      std::numeric_limits<double>::infinity();  // Time limit in seconds
  double tolerance = 1e-8;                      // Numerical tolerance
  bool useSlacks = false;    // Use slack variables in feasibility
  bool filterChecks = true;  // Use filtering for fast validity checks
  Exploration exploration = Exploration::ITERATIVE;  // Exploration strategy
  Branching branching = Branching::BINARY;           // Branching mode
  SplitScoring splitScoring = SplitScoring::MINMAX;  // Split scoring strategy
};

/**
 * @brief Dynamic programming solver for minimum-depth Voronoi decision trees.
 *
 * This class implements a dynamic programming algorithm that computes the
 * minimum number of branching operations needed in the worst case to locate
 * any input cost vector on a face of the Voronoi diagram. The result is used
 * to synthesize an optimal decision tree.
 */
class Dynprog {
public:
  /**
   * @brief Construct the dynamic programming solver.
   * @param voronoi The Voronoi diagram (must be built)
   * @param domain Optional constraints restricting the input domain
   */
  Dynprog(const Voronoi& voronoi, const Domain& domain = Domain());

  /**
   * @brief Run the dynamic programming algorithm.
   * @param params Algorithm parameters
   */
  void run(const DynprogParams& params = DynprogParams());

  /// Get algorithm parameters
  const DynprogParams& params() const { return params_; }

  /// Get algorithm statistics
  const DynprogStats& stats() const { return stats_; }

  /// Get algorithm logs
  const DynprogLogs& logs() const { return logs_; }

  /// Get the termination status
  DynprogStatus status() const { return status_; }

  /// Get a state by its index
  const State& state(Index id) const { return states_.at(id); }

  /// Get the root state index
  Index rootId() const { return rootId_; }

  const std::vector<Relation>& branchDirections() const {
    return branchDirections_;
  }

private:
  const Voronoi& voronoi_;  // Voronoi diagram reference
  Domain domain_;           // Input domain constraints
  DynprogParams params_;    // Algorithm parameters
  DynprogStats stats_;      // Algorithm statistics
  DynprogLogs logs_ = {};   // Progress logs

  DynprogStatus status_ = DynprogStatus::INVALID;  // Termination status
  std::vector<State> states_ = {};                 // All created states
  std::unordered_map<Cone, Index, ConeHasher> regionToStateId_ =
      {};                         // State lookup by region
  Index rootId_ = INVALID_INDEX;  // Root state index
  std::unique_ptr<Feasibility> feasibility_ = nullptr;  // Feasibility checker
  std::vector<Relation> branchDirections_ = {};  // Cached branching directions
  std::vector<std::vector<Relation>>
      positions_;                // Face/split position cache (LT/EQ/GT/RF)
  bool domainContainsOrigin_;    // Whether cost domain includes origin
  std::mt19937 rng_;             // Random generator for split scoring
  Clock::time_point startTime_;  // Algorithm start time
  Clock::time_point checkTime_;  // Last logging time

  void initRootState();
  void initPositions();

  State createState(State& state, Split& split, Relation cutDir);
  void buildState(Index stateId);
  void evaluateState(Index stateId, Index kSplits);
  void evaluateLowerBound(State& state);
  void evaluateScore(State& state, Split& split,
                     const std::vector<Index>& childFaceCounts);
  void evaluateValidity(State& state, Split& split,
                        bool externalStateCuts = false);
  bool checkChildFace(State& state, Split& split, Relation cutDir, Index faceId,
                      bool externalChildCuts = false);
  void updateStatus();

  std::optional<bool> inferValidity(State& state, Split& split);
  std::optional<bool> inferChildFace(State& state, Split& split,
                                     Relation cutDir, Index faceId);
  bool childContainsFaceCenter(State& state, Split& split, Relation cutDir,
                               Index faceId);

  std::vector<Index> getIterationRange() const;
  bool checkTimeLimit() const;
  void logHeader() const;
  void logProgress(const std::string& message = "");
  void logSave();
  void logFooter() const;
};

}  // namespace treeco

#endif  // TREECO_DYNPROG_HPP
