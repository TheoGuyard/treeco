#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <treeco/CLI.hpp>
#include <treeco/IO.hpp>
#include <treeco/LDTree.hpp>
#include <treeco/Types.hpp>
#include <vector>

using namespace treeco;

int main(int argc, char* argv[]) {
  Args args = parseArgs(argc, argv);

  std::ostream& out = std::cout;

  if (args.showHelp) {
    printUsage(argv[0], &out);
    return 0;
  }

  try {
    // ====================================================================
    // Load input data
    // ====================================================================

    out << "========================================\n";
    out << "Loading input data\n";
    out << "========================================\n";
    out << "\n";

    out << "Loading points file: " << args.pointsFile << "\n";
    std::vector<BinaryVector> points = readPoints(args.pointsFile);
    out << "  dim points: " << points[0].size() << "\n";
    out << "  num points: " << points.size() << "\n";

    Domain domain;
    if (!args.domainFile.empty()) {
      out << "Loading domain file: " << args.domainFile << "\n";
      domain = readDomain(args.domainFile);
      out << "  dim constr: " << std::to_string(std::get<0>(domain[0]).size())
          << "\n";
      out << "  num constr: " << domain.size() << "\n";

      if (points[0].size() != std::get<0>(domain[0]).size()) {
        std::string msg =
            "Dimension mismatch between points and domain constraints.";
        msg += "  dim points: " + std::to_string(points[0].size()) + "\n";
        msg +=
            "  dim constr: " + std::to_string(std::get<0>(domain[0]).size()) +
            "\n";
        throw std::runtime_error(msg);
      }
    }

    std::vector<RealVector> queries;
    if (!args.queriesFile.empty()) {
      out << "Loading queries file: " << args.queriesFile << "\n";
      queries = readQueries(args.queriesFile);
      out << "  dim queries: " << queries[0].size() << "\n";
      out << "  num queries: " << queries.size() << "\n";

      if (points[0].size() != queries[0].size()) {
        std::string msg = "Dimension mismatch between points and queries.";
        msg += "  dim points: " + std::to_string(points[0].size()) + "\n";
        msg += "  dim queries: " + std::to_string(queries[0].size()) + "\n";
        throw std::runtime_error(msg);
      }
    }

    // ====================================================================
    // Build LDT policy
    // ====================================================================

    out << "\n";
    out << "========================================\n";
    out << "Building LDT policy\n";
    out << "========================================\n";
    out << "\n";

    LDTree ldtree(points, domain);

    ldtree.build(args.verbose, &out, args.logInterval, args.logSave,
                 args.timeLimit, args.tolerance, args.deduplicate,
                 args.filterChecks, args.exploration, args.branching,
                 args.splitScoring);

    // ====================================================================
    // Perform queries
    // ====================================================================

    if (!queries.empty()) {
      out << "\n";
      out << "========================================\n";
      out << "Performing queries\n";
      out << "========================================\n";
      out << "\n";

      for (size_t q = 0; q < queries.size(); ++q) {
        const RealVector& cost = queries[q];
        std::vector<BinaryVector> solutions = ldtree.query(cost);

        // Compute objective values
        std::vector<double> values;
        values.reserve(solutions.size());
        for (const auto& sol : solutions) { values.push_back(dot(cost, sol)); }

        // Output results
        out << "query " << (q + 1) << "\n";
        out << "  cost : ";
        printVector(cost, &out);
        out << "\n";
        out << "  nsol : " << solutions.size() << "\n";
        out << "  sols : {";
        for (size_t s = 0; s < solutions.size(); ++s) {
          printVector(solutions[s], &out);
          if (s < solutions.size() - 1) out << ",";
        }
        out << "}\n";
        out << "  vals : {";
        for (size_t v = 0; v < values.size(); ++v) {
          out << values[v];
          if (v < values.size() - 1) out << ",";
        }
        out << "}\n";
        out << "\n";
      }
    } else {
      out << "\nNo queries to perform.\n";
    }

  } catch (const std::exception& e) {
    out << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
