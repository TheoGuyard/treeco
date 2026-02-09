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
#include <map>
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

/// Maximum possible or invalid depth value
constexpr Index MAX_DEPTH = std::numeric_limits<Index>::max();

/// Safely add positive integers without overflow
inline Index safeAdd(Index d, Index increment = 1) { return (d >= MAX_DEPTH - increment) ? MAX_DEPTH : d + increment; }

// ----------------------------------------------------------------------------
// Dynamic programming strategies
// ----------------------------------------------------------------------------

/// Exploration strategy for the dynamic programming algorithm
enum class Exploration {
  GREEDY,     // Single iteration with k=1 (potentially suboptimal, fast)
  ITERATIVE,  // Iterative scheme with k=1,2,...,kmax (optimal, early stopping)
  EXHAUSTIVE  // Single iteration with k=kmax (optimal, no early stopping)
};

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

inline Exploration stringToExplorationType(const std::string& str) {
  if (str == "greedy") return Exploration::GREEDY;
  if (str == "iterative") return Exploration::ITERATIVE;
  if (str == "exhaustive") return Exploration::EXHAUSTIVE;
  throw std::invalid_argument("Invalid exploration type: " + str +
                              ". Available options are: "
                              "greedy, iterative, exhaustive.");
}

/// Branching mode for decision tree nodes
enum class Branching {
  TERNARY,  // {<, =, >}
  BINARY    // {<, >} with tie-breaking for {=}
};

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

inline Branching stringToBranchingType(const std::string& str) {
  if (str == "binary") return Branching::BINARY;
  if (str == "ternary") return Branching::TERNARY;
  throw std::invalid_argument("Invalid branching type: " + str +
                              ". Available options are: "
                              "binary, ternary.");
}

/// Split scoring strategy for ordering candidate splits
enum class SplitScoring {
  VARIANCE,  // Variance-based deviation from equal split
  ENTROPY,   // Information gain (entropy reduction)
  MINMAX,    // Minimize maximum number of child faces
  NONE,      // All splits have equal score
  RANDOM     // Random scoring
};

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

// ----------------------------------------------------------------------------
// Dynamic programming sub-structures
// ----------------------------------------------------------------------------

/// Status of the dynamic programming algorithm
enum class DynprogStatus {
  INVALID,     // Not started, in progress, or terminated without valid tree
  SUBOPTIMAL,  // Found valid tree but optimality not certified
  OPTIMAL      // Found and certified optimal tree depth
};

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

/// Statistics from the dynamic programming algorithm
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

std::ostream& operator<<(std::ostream& oss, const DynprogStats& stats);

/// Logs from the dynamic programming algorithm
using DynprogLogs = std::vector<DynprogStats>;

/// Parameters for the dynamic programming algorithm
struct DynprogParams {
  bool verbose = false;                                        // Enable verbose logging
  std::ostream* outputStream = &std::cout;                     // Output stream for logs
  double logInterval = 5.0;                                    // Logging interval in seconds
  bool logSave = false;                                        // Save logs at each iteration
  double timeLimit = std::numeric_limits<double>::infinity();  // Time limit in seconds
  double tolerance = 1e-8;                                     // Numerical tolerance
  bool useSlacks = true;                                       // Use slack variables in feasibility checks
  bool strongChecks = true;                                    // Use strong checks while building states
  bool filterChecks = true;                                    // Use filtering for fast validity checks
  Exploration exploration = Exploration::ITERATIVE;            // Exploration strategy
  Branching branching = Branching::BINARY;                     // Branching mode
  SplitScoring splitScoring = SplitScoring::MINMAX;            // Split scoring strategy
};

// ----------------------------------------------------------------------------
// Dynamic programming state
// ----------------------------------------------------------------------------

/// Dynamic programming state split information
struct Split {
  Index id = INVALID_INDEX;                                // Split index in Voronoi diagram
  std::optional<bool> valid = std::nullopt;                // Whether the split is valid for this state
  double score = std::numeric_limits<double>::infinity();  // Split quality score
  std::array<Index, 3> childIds = {};                      // Child state indices
  bool isClosed = false;                                   // Whether the split is fully explored
};

/// Dynamic programming state
struct State {
  Cone region;                      // Polyhedral region of the state
  std::vector<Index> faceIds = {};  // Voronoi face indices in this region
  std::vector<Split> splits = {};   // Valid splits for this state
  Index lbHeight = 0;               // Lower bound on optimal subtree height
  Index ubHeight = MAX_DEPTH;       // Upper bound on optimal subtree height
  Index splitId = INVALID_INDEX;    // Best split index found
  bool isBuilt = false;             // Whether all splits have been explored
  bool isClosed = false;            // Whether state is closed (leaf, pruned, or optimal)

  /// Get the depth of this state in the tree
  Index depth() const { return static_cast<Index>(region.numCuts()); }

  /// Check if this state is a leaf (single face or no valid splits)
  bool isLeaf() const { return faceIds.size() <= 1 || splits.size() == 0; }
};

// ----------------------------------------------------------------------------
// Dynamic programming method
// ----------------------------------------------------------------------------

/// Dynamic programming memoization structure
using DynprogMemo = std::unordered_map<Cone, Index, ConeHasher>;

/// Dynamic programming solver
class Dynprog {
public:
  /// Construct the dynamic programming solver
  Dynprog(const Voronoi& voronoi, const Domain& domain = Domain());

  /// Run the dynamic programming algorithm
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

  /// Get branching directions
  const std::vector<Relation>& branchDirections() const { return branchDirections_; }

private:
  const Voronoi& voronoi_;  // Voronoi diagram reference
  Domain domain_;           // Input domain constraints
  DynprogParams params_;    // Algorithm parameters
  DynprogStats stats_;      // Algorithm statistics
  DynprogLogs logs_ = {};   // Progress logs

  DynprogStatus status_ = DynprogStatus::INVALID;       // Termination status
  std::vector<State> states_ = {};                      // All created states
  DynprogMemo regionToStateId_ = {};                    // State lookup by region
  Index rootId_ = INVALID_INDEX;                        // Root state index
  std::unique_ptr<Feasibility> feasibility_ = nullptr;  // Feasibility checker
  std::vector<Relation> branchDirections_ = {};         // Cached branching directions
  std::vector<std::vector<Relation>> positions_;        // Face/split position cache (LT/EQ/GT/RF)
  Clock::time_point startTime_;                         // Algorithm start time
  Clock::time_point checkTime_;                         // Last logging time

  void initRootState();
  void initPositions();

  Index getOrCreateChildId(Index parentId, Index splitPos, Relation cutDir);
  State createState(const State& parent, const Split& split, Relation cutDir);
  State dummyChild(const State& parent, const Split& split, Relation cutDir, const Cone& region);
  void buildState(Index stateId);
  void evaluateState(Index stateId, Index kSplits);
  void evaluateLowerBound(State& state);
  void evaluateScore(State& state, Split& split, const std::vector<Index>& childFaceCounts);
  bool checkChildFace(const State& state, const Split& split, Relation cutDir, Index faceId,
                      bool externalChildCuts = false);
  void updateStatus();

  std::optional<bool> inferChildFace(const State& state, const Split& split, Relation cutDir, Index faceId);
  bool childContainsFaceCenter(const State& state, const Split& split, Relation cutDir, Index faceId);

  std::vector<Index> getIterationRange() const;
  bool checkTimeLimit() const;
  void logHeader() const;
  void logProgress(const std::string& message = "");
  void logSave();
  void logFooter() const;
};

}  // namespace treeco

#endif  // TREECO_DYNPROG_HPP
