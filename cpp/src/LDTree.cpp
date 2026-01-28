#include "treeco/LDTree.hpp"

namespace treeco {

LDTree::LDTree(const std::string &filePoints, const std::string &fileDomain)
    : LDTree(readPoints(filePoints),
             fileDomain.empty() ? Domain() : readDomain(fileDomain)) {}

LDTree::LDTree(const std::vector<BinaryVector> &points, const Domain &domain)
    : domain_(domain), voronoi_(scaleBinarySet(points)), tree_(voronoi_) {}

void LDTree::build(bool verbose, std::ostream *outputStream, double logInterval,
                   double timeLimit, double tolerance, bool deduplicate,
                   bool filterChecks, Exploration exploration,
                   Branching branching, LowerBounding lowerBounding,
                   Positioning positioning, SplitSelection splitSelection,
                   SplitScoring splitScoring, Index randomSeed) {

  auto startTime = Clock::now();

  if (verbose) {
    *outputStream << "Building LDTree...\n";
  }

  // Step 1: Build the Voronoi diagram structure
  VoronoiParams voronoiParams = VoronoiParams{
      verbose,   outputStream, logInterval, timeLimit - elapsedTime(startTime),
      tolerance, deduplicate};
  voronoi_.build(voronoiParams);

  // Step 2: Run dynamic programming to find the minimum tree depth
  DynprogParams dynprogParams =
      DynprogParams{verbose,        outputStream,
                    logInterval,    timeLimit - elapsedTime(startTime),
                    tolerance,      filterChecks,
                    exploration,    branching,
                    lowerBounding,  positioning,
                    splitSelection, splitScoring,
                    randomSeed};
  Dynprog dynprog = Dynprog(voronoi_, domain_);
  dynprog.run(dynprogParams);

  // Step 3: Synthetize tree structure from the dynamic programming results
  TreeParams treeParams =
      TreeParams{verbose, outputStream, logInterval,
                 timeLimit - elapsedTime(startTime), tolerance};
  tree_.synthetize(dynprog, treeParams);

  // Record statistics
  stats_.buildTime = elapsedTime(startTime);

  if (verbose) {
    std::ostream &out = *(outputStream);
    out << "LDTree constructed\n";
    out << "  total time  : " << std::fixed << std::setprecision(4)
        << stats_.buildTime << "\n";
    out << "  tree size   : " << tree_.size() << "\n";
    out << "  tree width  : " << tree_.width() << "\n";
    out << "  tree depth  : " << tree_.depth() << "\n";
  }
};

std::vector<BinaryVector> LDTree::query(const RealVector &cost,
                                        bool checkDomain) const {
  if (checkDomain) {
    for (const auto &[a, b, rel] : domain_) {
      double value = dot(a, cost) + b;
      if ((rel == Relation::LT && value >= 0.0) ||
          (rel == Relation::LE && value > 0.0) ||
          (rel == Relation::EQ && std::abs(value) > 0.0) ||
          (rel == Relation::GE && value < 0.0) ||
          (rel == Relation::GT && value <= 0.0) || (rel == Relation::RF)) {
        throw std::runtime_error("Cost vector is outside the defined domain.");
      }
    }
  }
  return unscaleBinarySet(tree_.query(cost));
}

void LDTree::pprint(bool tightDisplay, std::ostream *outputStream) const {
  return tree_.pprint(tightDisplay, outputStream);
}

void LDTree::flatten(const std::string filepath, const std::string doc,
                     bool benchmarkMode) const {

  Index dimPoints = voronoi_.dimPoints();

  std::ostringstream out;

  // Header
  out << "#include <stdio.h>\n";
  out << "#include <stdlib.h>\n";
  if (benchmarkMode) {
    out << "#include <time.h>\n";
  }
  out << "\n";
  out << "#define DIM " << dimPoints << "\n";
  out << "\n";
  out << "\n";

  // Query function
  out << "static inline const int *query(const double *c) {\n";
  generateNodeCode(0, out, 1);
  out << "}\n";
  out << "\n";

  // Main function
  if (benchmarkMode) {
    generateBenchmarkCode(out, doc);
  } else {
    generateNormalCode(out, doc);
  }
  out << "\n";

  // Save code to file
  std::ofstream file(filepath);
  file << out.str();
}

void LDTree::generateNormalCode(std::ostringstream &out,
                                const std::string &doc) const {
  out << "int main(int argc, char **argv) {\n";
  out << "    if (argc != DIM + 1) {\n";
  out << "        fprintf(stderr, \"Usage: %s " << doc << " \\n\", argv[0]);\n";
  out << "        return 1;\n";
  out << "    }\n";
  out << "\n";
  out << "    double c[DIM];\n";
  out << "    for (int i = 0; i < DIM; i++) {\n";
  out << "        c[i] = atof(argv[i + 1]);\n";
  out << "    }\n";
  out << "\n";
  out << "    const int *sol = query(c);\n";
  out << "\n";
  out << "    printf(\"Solution: [\");\n";
  out << "    for (int i = 0; i < DIM; i++) {\n";
  out << "        printf(\"%d\", sol[i]);\n";
  out << "        if (i + 1 < DIM) printf(\", \");\n";
  out << "    }\n";
  out << "    printf(\"]\\n\");\n";
  out << "    return 0;\n";
  out << "}\n";
}

void LDTree::generateBenchmarkCode(std::ostringstream &out,
                                   const std::string &doc) const {
  out << "int main(int argc, char **argv) {\n";
  out << "    if (argc != 3) {\n";
  out << "        fprintf(stderr, \"Usage: %s <runs> <seed> " << doc
      << " \\n\", argv[0]);\n";
  out << "        return 1;\n";
  out << "    }\n";
  out << "\n";
  out << "    int runs = atoi(argv[1]);\n";
  out << "    int seed = atoi(argv[2]);\n";
  out << "\n";
  out << "    srand(seed);\n";
  out << "    double **costs = malloc(runs * sizeof(double*));\n";
  out << "    for (int i = 0; i < runs; ++i) {\n";
  out << "        costs[i] = malloc(DIM * sizeof(double));\n";
  out << "        for (int j = 0; j < DIM; ++j) {\n";
  out << "            costs[i][j] = (rand() + 1.0) / RAND_MAX;\n";
  out << "        }\n";
  out << "    }\n";
  out << "\n";
  out << "    struct timespec start, end;\n";
  out << "    clock_gettime(CLOCK_MONOTONIC, &start);\n";
  out << "    for (int i = 0; i < runs; ++i) {\n";
  out << "        const int *out = query(costs[i]);\n";
  out << "    }\n";
  out << "    clock_gettime(CLOCK_MONOTONIC, &end);\n";
  out << "\n";
  out << "    for (int i = 0; i < runs; ++i) {\n";
  out << "        free(costs[i]);\n";
  out << "    }\n";
  out << "    free(costs);\n";
  out << "\n";
  out << "    long double elapsed = (end.tv_sec - start.tv_sec) +\n";
  out << "                     (long double)(end.tv_nsec - start.tv_nsec) / "
         "1e9;\n";
  out << "    printf(\"  runs       : %d\\n\", runs);\n";
  out << "    printf(\"  seed       : %d\\n\", seed);\n";
  out << "    printf(\"  total_time : %.12Lf\\n\", elapsed);\n";
  out << "    printf(\"  avg_time   : %.12Lf\\n\", elapsed / runs);\n";
  out << "    return 0;\n";
  out << "}\n";
}

void LDTree::generateNodeCode(Index nodeId, std::ostringstream &out,
                              int indent) const {
  const Node &node = tree_.node(nodeId);
  std::string pad(indent * 4, ' ');
  if (node.type == NodeType::LEAF) {
    BinaryVector point = unscaleBinary(voronoi_.point(node.pointsIds[0]));
    out << pad << "static int sol[DIM] = {";
    for (Index i = 0; i < point.size(); i++) {
      out << point[i];
      out << ((i + 1 < point.size()) ? ", " : "};\n");
    }
    out << pad << "return sol;\n";
    return;
  }
  out << pad << "if (";
  generateDotCode(out, voronoi_.split(node.splitId));
  out << " <= 0.0) {\n";
  generateNodeCode(node.childIds.at(Relation::LT), out, indent + 1);
  out << pad << "} else {\n";
  generateNodeCode(node.childIds.at(Relation::GT), out, indent + 1);
  out << pad << "}\n";
}

void LDTree::generateDotCode(std::ostringstream &out,
                             const TernaryVector &split) const {
  bool first = true;
  for (Index i = 0; i < split.size(); i++) {
    if (split[i] == -1) {
      if (!first) {
        out << " - ";
      }
      out << "c[" << i << "]";
      first = false;
    } else if (split[i] == 0) {
      continue;
    } else if (split[i] == 1) {
      if (!first) {
        out << " + ";
      }
      out << "c[" << i << "]";
      first = false;
    } else {
      throw std::runtime_error("Invalid entry value.");
    }
  }
}

std::ostream &operator<<(std::ostream &oss, const LDTree &ldtree) {
  oss << "LDTree\n";
  oss << "  dim points: " << ldtree.voronoi().dimPoints() << "\n";
  oss << "  num points: " << ldtree.voronoi().numPoints() << "\n";
  oss << "  max depth : " << ldtree.tree().depth();
  return oss;
}

} // namespace treeco