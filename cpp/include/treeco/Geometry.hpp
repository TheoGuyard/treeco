/**
 * @file Geometry.hpp
 * @brief Geometric structures and utilities for hyperplane arrangements.
 *
 * This header defines cuts, cones, and geometric utility functions used
 * for representing and manipulating polyhedral regions in the decision tree.
 */

#ifndef TREECO_GEOMETRY_HPP
#define TREECO_GEOMETRY_HPP

#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "treeco/Types.hpp"

namespace treeco {

// ----------------------------------------------------------------------------
// Relation types
// ----------------------------------------------------------------------------

/**
 * @brief Relation with respect to a hyperplane defined by <a,x> + b.
 *
 * Represents the position of a point relative to a hyperplane, used for
 * defining constraints and branching decisions in the decision tree.
 */
enum class Relation : int8_t {
  LT = -2,  // <a,c> + b < 0 (strict less than)
  LE = -1,  // <a,c> + b <= 0 (less than or equal)
  EQ = 0,   // <a,c> + b = 0 (equal)
  GE = 1,   // <a,c> + b >= 0 (greater than or equal)
  GT = 2,   // <a,c> + b > 0 (strict greater than)
  RT = 10,  // <a,c> + b <= +inf (always true, unconstrained)
  RF = 20,  // <a,c> + b >= +inf (always false, infeasible)
};

/**
 * @brief Convert a Relation enum to its string representation.
 * @param type The relation type
 * @return String representation of the relation
 */
inline std::string relationTypeToString(Relation type) {
  switch (type) {
    case Relation::LT:
      return "LT";
    case Relation::LE:
      return "LE";
    case Relation::EQ:
      return "EQ";
    case Relation::GE:
      return "GE";
    case Relation::GT:
      return "GT";
    case Relation::RT:
      return "RT";
    case Relation::RF:
      return "RF";
    default:
      return "Unknown";
  }
}

/**
 * @brief Parse a string to a Relation enum.
 * @param str The string representation (LT, LE, EQ, GE, GT, RT, RF)
 * @return The corresponding Relation enum value
 * @throws std::invalid_argument if the string is not a valid relation
 */
inline Relation stringToRelationType(const std::string& str) {
  if (str == "LT") return Relation::LT;
  if (str == "LE") return Relation::LE;
  if (str == "EQ") return Relation::EQ;
  if (str == "GE") return Relation::GE;
  if (str == "GT") return Relation::GT;
  if (str == "RT") return Relation::RT;
  if (str == "RF") return Relation::RF;
  throw std::invalid_argument(
      "Invalid relation type string: " + str +
      ". Available options are: " + "LT, LE, EQ, GE, GT, RT, RF.");
}

/// Lookup table for reducing pairs of relations
const std::map<std::set<Relation>, Relation> relationReductions = {
    {{Relation::LT, Relation::LE}, Relation::LT},
    {{Relation::LE, Relation::EQ}, Relation::EQ},
    {{Relation::LE, Relation::GE}, Relation::EQ},
    {{Relation::EQ, Relation::GE}, Relation::EQ},
    {{Relation::GE, Relation::GT}, Relation::GT},
};

/**
 * @brief Reduce a pair of relations to their intersection.
 * @param r1 First relation
 * @param r2 Second relation
 * @return The reduced relation representing the intersection
 */
inline Relation reduceRelationPair(Relation r1, Relation r2) {
  if (r1 == r2) {
    return r1;
  } else if (r1 == Relation::RT) {
    return r2;
  } else if (r2 == Relation::RT) {
    return r1;
  } else if (r1 == Relation::RF || r2 == Relation::RF) {
    return Relation::RF;
  } else {
    std::set<Relation> key = {r1, r2};
    auto it = relationReductions.find(key);
    return (it != relationReductions.end()) ? it->second : Relation::RF;
  }
}

/**
 * @brief Reduce a vector of relations to their intersection.
 * @param relations Vector of relations to reduce
 * @return The reduced relation (RT if empty, RF if inconsistent)
 */
inline Relation reduceRelations(const std::vector<Relation>& relations) {
  if (relations.empty()) return Relation::RT;
  return std::accumulate(std::next(relations.begin()), relations.end(),
                         relations.front(), reduceRelationPair);
}

// ----------------------------------------------------------------------------
// Constraint and Domain types
// ----------------------------------------------------------------------------

/// Tuple (a,b,rel) representing a linear constraint <a,x> + b rel 0
using ConstrData = std::tuple<RealVector, double, Relation>;

/// Collection of linear constraints defining a polyhedral domain
using Domain = std::vector<ConstrData>;

/// Positive orthant {x in R^n | x[i] > 0 for all i} of given dimension
inline Domain positiveOrthant(Index dim) {
  Domain domain;
  for (Index i = 0; i < dim; ++i) {
    RealVector a(dim, 0.0);
    a[i] = 1.0;
    domain.emplace_back(a, 0.0, Relation::GT);
  }
  return domain;
}

/// Non-negative orthant {x in R^n | x[i] >= 0 for all i} of given dimension
inline Domain nonNegativeOrthant(Index dim) {
  Domain domain;
  for (Index i = 0; i < dim; ++i) {
    RealVector a(dim, 0.0);
    a[i] = 1.0;
    domain.emplace_back(a, 0.0, Relation::GE);
  }
  return domain;
}

/// Negative orthant {x in R^n | x[i] < 0 for all i} of given dimension
inline Domain negativeOrthant(Index dim) {
  Domain domain;
  for (Index i = 0; i < dim; ++i) {
    RealVector a(dim, 0.0);
    a[i] = 1.0;
    domain.emplace_back(a, 0.0, Relation::LT);
  }
  return domain;
}

/// Non-positive orthant {x in R^n | x[i] <= 0 for all i} of given dimension
inline Domain nonPositiveOrthant(Index dim) {
  Domain domain;
  for (Index i = 0; i < dim; ++i) {
    RealVector a(dim, 0.0);
    a[i] = 1.0;
    domain.emplace_back(a, 0.0, Relation::LE);
  }
  return domain;
}

// ----------------------------------------------------------------------------
// Geometrical structures
// ----------------------------------------------------------------------------

/**
 * @brief A cut representing a half-space constraint <v,x> {<,<=,=,>=,>} 0.
 *
 * The cut is defined by the index of its associated normal vector v (from a
 * Voronoi diagram's split pool) and the relation type defining the half-space.
 */
struct Cut {
  Index hid = INVALID_INDEX;    // Index of the normal vector in the split pool
  Relation dir = Relation::RT;  // Direction/relation of the cut

  /**
   * @brief Construct a cut.
   * @param hid Index of the hyperplane normal vector
   * @param dir Relation type (LT, LE, EQ, GE, GT)
   */
  Cut(Index hid, Relation dir);

  /// Compute a hash value for this cut
  std::size_t hash() const;

  /// Get string representation of the cut
  std::string toString() const;

  /// Check equality with another cut
  bool operator==(const Cut& other) const;

  /// Check inequality with another cut
  bool operator!=(const Cut& other) const;

  /// Lexicographic ordering for use in ordered containers
  bool operator<(const Cut& other) const;
};

/// Hasher functor for Cut objects (for use in unordered containers)
struct CutHasher {
  std::size_t operator()(const Cut& c) const { return c.hash(); }
};

/**
 * @brief A polyhedral cone defined as the intersection of half-spaces.
 *
 * The cone is represented as a set of cuts, where each cut defines a
 * half-space boundary. The cone is the intersection of all these half-spaces.
 */
struct Cone {
public:
  /// Construct an empty cone (full space)
  Cone();

  /**
   * @brief Construct a cone from a set of cuts.
   * @param cuts The cuts defining the cone boundaries
   */
  Cone(const std::set<Cut>& cuts);

  /**
   * @brief Construct a cone from a vector of cuts.
   * @param cuts The cuts defining the cone boundaries
   */
  Cone(const std::vector<Cut>& cuts);

  /**
   * @brief Add a cut boundary to the cone.
   * @param cut The cut to add
   */
  void addCut(const Cut& cut);

  /**
   * @brief Add a cut boundary to the cone.
   * @param hid Index of the hyperplane normal vector
   * @param dir Relation type
   */
  void addCut(Index hid, Relation dir);

  /// Get the cone cuts
  const std::set<Cut>& cuts() const;

  /// Get the number of cut boundaries
  Index numCuts() const;

  /// Check if the cone is the full space (no boundaries)
  bool isFullSpace() const;

  /// Check whether the cone is open
  bool isOpen() const;

  /// Check whether the cone contains the origin
  bool containsOrigin() const;

  /// Remove all boundaries (reset to full space)
  void clear();

  /**
   * @brief Create a refined cone by adding a new cut.
   * @param cut The cut to add
   * @return A new cone with the additional cut
   */
  Cone refine(const Cut& cut) const;

  /**
   * @brief Get the interior of the cone (cuts set to strict inequalities).
   * @return A new cone with all cuts as strict inequalities
   */
  Cone interior() const;

  /// Compute a hash value for this cone (for memoization)
  std::size_t hash() const;

  /// Check equality with another cone
  bool operator==(const Cone& other) const;

private:
  /// Set of cuts defining the cone boundaries
  std::set<Cut> cuts_ = {};

  /// Whether the cone is open
  bool isOpen_ = true;

  /// Whether the cone contains the origin
  bool containsOrigin_ = true;
};

/// Hasher functor for Cone objects (for use in unordered containers)
struct ConeHasher {
  std::size_t operator()(const Cone& c) const { return c.hash(); }
};

/// Stream output operator for Cone objects
std::ostream& operator<<(std::ostream& oss, const Cone& cone);

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

/**
 * @brief Normalize a ternary vector so that its first non-zero entry is +1.
 * @param v The vector to normalize (modified in place)
 */
void normalize(TernaryVector& v);

/**
 * @brief Compute the dot product of a ternary vector and a simplex vector.
 * @param p Ternary vector in {-1,0,+1}^n
 * @param v Simplex vector in {-1,+1}^n
 * @return Integer dot product
 */
int dot(const TernaryVector& p, const SimplexVector& v);

/**
 * @brief Compute the dot product of a ternary vector and a real vector.
 * @param p Ternary vector in {-1,0,+1}^n
 * @param v Real vector in R^n
 * @return Real-valued dot product
 */
double dot(const TernaryVector& p, const RealVector& v);

/**
 * @brief Compute the dot product of a real vector and a binary vector.
 * @param p Real vector in R^n
 * @param v Binary vector in {0,1}^n
 * @return Real-valued dot product (objective value)
 */
double dot(const RealVector& p, const BinaryVector& v);

/**
 * @brief Compute the dot product of real vectors.
 * @param p Real vector in R^n
 * @param v Real vector in R^n
 * @return Real-valued dot product (objective value)
 */
double dot(const RealVector& p, const RealVector& v);

/**
 * @brief Compute the bisector hyperplane between two simplex points.
 * @param p1 First point in {-1,+1}^n
 * @param p2 Second point in {-1,+1}^n
 * @return Ternary vector representing the separating hyperplane normal
 */
TernaryVector bisector(const SimplexVector& p1, const SimplexVector& p2);

/**
 * @brief Map a binary vector to simplex coordinates: x -> 2x - 1.
 * @param x Binary vector in {0,1}^n
 * @return Simplex vector in {-1,+1}^n
 */
SimplexVector scaleBinary(const BinaryVector& x);

/**
 * @brief Map a simplex vector to binary coordinates: x -> (x + 1) / 2.
 * @param x Simplex vector in {-1,+1}^n
 * @return Binary vector in {0,1}^n
 */
SimplexVector unscaleBinary(const SimplexVector& x);

/**
 * @brief Map a set of binary vectors to simplex coordinates.
 * @param X Set of binary vectors in {0,1}^n
 * @return Set of simplex vectors in {-1,+1}^n
 */
std::vector<SimplexVector> scaleBinarySet(const std::vector<BinaryVector>& X);

/**
 * @brief Map a set of simplex vectors to binary coordinates.
 * @param X Set of simplex vectors in {-1,+1}^n
 * @return Set of binary vectors in {0,1}^n
 */
std::vector<SimplexVector> unscaleBinarySet(
    const std::vector<SimplexVector>& X);

}  // namespace treeco

#endif  // TREECO_GEOMETRY_HPP
