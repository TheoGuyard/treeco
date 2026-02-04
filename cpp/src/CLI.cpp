#include <treeco/CLI.hpp>

namespace treeco {

Args parseArgs(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "--help" || arg == "-h") {
      args.showHelp = true;
    } else if (arg == "--points" && i + 1 < argc) {
      args.pointsFile = argv[++i];
    } else if (arg == "--domain" && i + 1 < argc) {
      args.domainFile = argv[++i];
    } else if (arg == "--queries" && i + 1 < argc) {
      args.queriesFile = argv[++i];
    } else if (arg == "--verbose" || arg == "-v") {
      args.verbose = true;
    } else if (arg == "--log-interval" && i + 1 < argc) {
      args.logInterval = std::stod(argv[++i]);
    } else if (arg == "--log-save") {
      args.logSave = true;
    } else if (arg == "--time-limit" && i + 1 < argc) {
      args.timeLimit = std::stod(argv[++i]);
    } else if (arg == "--tolerance" && i + 1 < argc) {
      args.tolerance = std::stod(argv[++i]);
    } else if (arg == "--no-deduplicate") {
      args.deduplicate = false;
    } else if (arg == "--use-slacks") {
      args.useSlacks = true;
    } else if (arg == "--no-filter-checks") {
      args.filterChecks = false;
    } else if (arg == "--exploration" && i + 1 < argc) {
      args.exploration = stringToExplorationType(argv[++i]);
    } else if (arg == "--branching" && i + 1 < argc) {
      args.branching = stringToBranchingType(argv[++i]);
    } else if (arg == "--split-scoring" && i + 1 < argc) {
      args.splitScoring = stringToSplitScoringType(argv[++i]);
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      std::cerr << "Use --help for usage information.\n";
      std::exit(1);
    }
  }

  // Check that required arguments are provided
  if (args.pointsFile.empty() && !args.showHelp) {
    std::cerr << "Error: --points argument is required.\n";
    std::cerr << "Use --help for usage information.\n";
    std::exit(1);
  }

  return args;
}

void printUsage(const char* progName, std::ostream* outputStream) {
  std::ostream& out = (outputStream != nullptr) ? *outputStream : std::cout;
  out << "Usage: " << progName << " [options]\n\n";
  out << "Input files:\n";
  out << "  --points <file>       Points file (required, default: none)\n";
  out << "  --domain <file>       Domain file (optional, default: none)\n";
  out << "  --queries <file>      Queries file (optional, default: none)\n";
  out << "\n";
  out << "General parameters:\n";
  out << "  --verbose             Enable verbose output\n";
  out << "  --log-interval <sec>  Log interval in seconds (default: 5.0)\n";
  out << "  --log-save            Save logs at each iteration (default: "
         "true)\n";
  out << "  --time-limit <sec>    Time limit in seconds (default: infinity)\n";
  out << "  --tolerance <value>   Tolerance for zero (default: 1e-8)\n";
  out << "\n";
  out << "Voronoi parameters:\n";
  out << "  --no-deduplicate      Disable point deduplication\n";
  out << "\n";
  out << "Dynamic programming parameters:\n";
  out << "  --use-slacks              Use of slack variables in feasibility "
         "checks\n";
  out << "  --no-filter-checks        Disable filter checks optimization\n";
  out << "  --exploration <mode>      Exploration: greedy, iterative, "
         "exhaustive (default: iterative)\n";
  out << "  --branching <mode>        Branching: binary, ternary (default: "
         "binary)\n";
  out << "  --lower-bounding <mode>   Lower bounding: fixed, backtrack "
         "(default: backtrack)\n";
  out << "  --split-scoring <mode>    Split scoring: variance, entropy, "
         "minmax, balance, none, random (default: variance)\n";
  out << "\n";
  out << "Other:\n";
  out << "  --help                Show this help message\n";
}

}  // namespace treeco
