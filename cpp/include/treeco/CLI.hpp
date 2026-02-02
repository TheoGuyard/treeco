/**
 * @file CLI.hpp
 * @brief Command-line interface utilities.
 *
 * This header provides argument parsing and help display for the treeco
 * command-line executable.
 */

#ifndef TREECO_CLI_HPP
#define TREECO_CLI_HPP

#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <treeco/Dynprog.hpp>
#include <vector>

namespace treeco {

/**
 * @brief Command-line arguments structure.
 */
struct Args {
  std::string pointsFile = "";   // Path to points file
  std::string domainFile = "";   // Path to domain file (optional)
  std::string queriesFile = "";  // Path to queries file (optional)

  bool verbose = true;       // Enable verbose output
  double logInterval = 5.0;  // Logging interval in seconds
  bool logSave = true;       // Save logs at each iteration
  double timeLimit = std::numeric_limits<double>::infinity();  // Time limit
  double tolerance = 1e-8;  // Numerical tolerance

  bool deduplicate = true;  // Remove duplicate points

  bool filterChecks = true;  // Use filtering for validity checks
  Exploration exploration = Exploration::ITERATIVE;  // Exploration strategy
  Branching branching = Branching::BINARY;           // Branching mode
  SplitScoring splitScoring = SplitScoring::MINMAX;  // Split scoring strategy

  bool showHelp = false;  // Display help message
};

/**
 * @brief Parse command-line arguments.
 * @param argc Argument count
 * @param argv Argument values
 * @return Parsed arguments structure
 */
Args parseArgs(int argc, char* argv[]);

/**
 * @brief Print usage information.
 * @param progName Program name (argv[0])
 * @param outputStream Output stream for the help text
 */
void printUsage(const char* progName, std::ostream* outputStream);

}  // namespace treeco

#endif  // TREECO_CLI_HPP