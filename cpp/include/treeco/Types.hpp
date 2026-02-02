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
template <typename T>
inline void printVector(const std::vector<T>& vec,
                        std::ostream* out = &std::cout) {
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

}  // namespace treeco

#endif  // TREECO_TYPES_HPP
