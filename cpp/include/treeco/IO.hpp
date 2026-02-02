/**
 * @file IO.hpp
 * @brief Input/output utilities for reading and writing data files.
 *
 * This header provides functions for reading and writing points, domain
 * constraints, and query vectors from/to text files.
 */

#ifndef TREECO_IO_HPP
#define TREECO_IO_HPP

#include <fstream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "treeco/Dynprog.hpp"
#include "treeco/Geometry.hpp"
#include "treeco/Types.hpp"

namespace treeco {

/**
 * @brief Read feasible points from a file.
 * @param filepath Path to the points file
 * @return Vector of binary vectors representing feasible solutions
 */
std::vector<BinaryVector> readPoints(const std::string& filepath);

/**
 * @brief Write feasible points to a file.
 * @param filepath Output file path
 * @param points Vector of binary vectors to write
 */
void writePoints(const std::string& filepath,
                 const std::vector<BinaryVector>& points);

/**
 * @brief Read domain constraints from a file.
 * @param filepath Path to the domain file
 * @return Vector of constraint tuples (coefficients, relation, rhs)
 */
Domain readDomain(const std::string& filepath);

/**
 * @brief Write domain constraints to a file.
 * @param filepath Output file path
 * @param domain Vector of constraint tuples to write
 */
void writeDomain(const std::string& filepath, const Domain& domain);

/**
 * @brief Read query cost vectors from a file.
 * @param filepath Path to the queries file
 * @return Vector of real-valued cost vectors
 */
std::vector<RealVector> readQueries(const std::string& filepath);

/**
 * @brief Write query cost vectors to a file.
 * @param filepath Output file path
 * @param queries Vector of cost vectors to write
 */
void writeQueries(const std::string& filepath,
                  const std::vector<RealVector>& queries);

}  // namespace treeco

#endif  // TREECO_IO_HPP