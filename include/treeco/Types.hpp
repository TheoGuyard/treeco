/**
 * @file Types.hpp
 * @brief Core type definitions and utilities for the treeco library.
 * 
 * This header defines the fundamental types used throughout the library,
 * including vector types, index types, relation enums, and utility functions.
 */

#ifndef TREECO_TYPES_HPP
#define TREECO_TYPES_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace treeco {

/// Clock type for timing operations
using Clock = std::chrono::high_resolution_clock; 

/**
 * @brief Compute elapsed time since a start time point.
 * @param start The starting time point
 * @return Elapsed time in seconds
 */
inline double elapsedTime(const std::chrono::time_point<Clock>& start) {
    std::chrono::time_point<Clock> end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

/// Index type for array indexing
using Index = std::size_t;

/// Sentinel value for invalid indices
constexpr Index INVALID_INDEX = std::numeric_limits<Index>::max();

using BinaryVector = std::vector<int8_t>;   // Binary vector in {0,1}^n
using TernaryVector = std::vector<int8_t>;  // Ternary vector in {-1,0,+1}^n
using SimplexVector = std::vector<int8_t>;  // Simplex vector in {-1,+1}^n
using RealVector = std::vector<double>;     // Real-valued vector in R^n  

/**
 * @brief Print a vector to an output stream.
 * @tparam T The element type of the vector
 * @param vec The vector to print
 * @param out Output stream (default: std::cout)
 */
template<typename T>
inline void printVector(const std::vector<T>& vec, std::ostream* out = &std::cout) {
    (*out) << "[";
    for (Index i = 0; i < vec.size(); ++i) {
        if constexpr (std::is_same_v<T, int8_t>) {
            (*out) << static_cast<int>(vec[i]);
        } else {
            (*out) << vec[i];
        }
        if (i < vec.size() - 1) (*out) << ",";
    }
    (*out) << "]";
}


/**
 * @brief Relation with respect to a hyperplane defined by <a,x> + b.
 * 
 * Represents the position of a point relative to a hyperplane, used for
 * defining constraints and branching decisions in the decision tree.
 */
enum class Relation : int8_t {
    LT = -2,    // <a,c> + b < 0 (strict less than)
    LE = -1,    // <a,c> + b <= 0 (less than or equal)
    EQ = 0,     // <a,c> + b = 0 (equal)
    GE = 1,     // <a,c> + b >= 0 (greater than or equal)
    GT = 2,     // <a,c> + b > 0 (strict greater than)
    RT = 10,    // <a,c> + b <= +inf (always true, unconstrained)
    RF = 20,    // <a,c> + b >= +inf (always false, infeasible)
};

/**
 * @brief Convert a Relation enum to its string representation.
 * @param type The relation type
 * @return String representation of the relation
 */
inline std::string relationTypeToString(Relation type) {
    switch (type) {
        case Relation::LT: return "LT";
        case Relation::LE: return "LE";
        case Relation::EQ: return "EQ";
        case Relation::GE: return "GE";
        case Relation::GT: return "GT";
        case Relation::RT: return "RT";
        case Relation::RF: return "RF";
        default: return "Unknown";
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
        "Invalid relation type string: " + str + ". Available options are: "+
        "LT, LE, EQ, GE, GT, RT, RF."
    );
}

/// Tuple (a,b,rel) representing a linear constraint <a,x> + b rel 0
using ConstrData = std::tuple<RealVector, double, Relation>;

/// Collection of linear constraints defining a polyhedral domain
using Domain = std::vector<ConstrData>;

} // namespace treeco

#endif // TREECO_TYPES_HPP
