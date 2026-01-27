#include "treeco/IO.hpp"

#include <iomanip>

namespace treeco {

// ============================================================================
// Points I/O
// ============================================================================

std::vector<BinaryVector> readPoints(const std::string &filepath) {

  std::ifstream infile(filepath);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file.");
  }

  std::string header;
  if (!std::getline(infile, header)) {
    throw std::runtime_error("File is empty.");
  }

  std::istringstream ss(header);
  std::string keyword;
  Index numPoints;
  Index dimPoints;

  if (!(ss >> keyword >> dimPoints >> numPoints)) {
    throw std::runtime_error("Invalid header format.");
  } else if (keyword != "points") {
    throw std::runtime_error("Invalid header format.");
  }

  if (numPoints <= 0 || dimPoints <= 0) {
    throw std::runtime_error(
        "Number and dimension of points must be positive.");
  }

  std::vector<BinaryVector> points;
  points.reserve(numPoints);

  std::string line;
  Index pointsCount = 0;
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;

    std::istringstream row_stream(line);
    BinaryVector point;
    point.reserve(dimPoints);

    int entry;
    while (row_stream >> entry) {
      if (entry != 0 && entry != 1) {
        throw std::runtime_error("Point entries must be 0 or 1.");
      }
      point.push_back(entry);
    }

    if (point.size() != dimPoints) {
      throw std::runtime_error("Incorrect point dimension.");
    }

    points.push_back(std::move(point));
    pointsCount++;
  }

  if (pointsCount != numPoints) {
    throw std::runtime_error("Incorrect number of points");
  }

  return points;
}

void writePoints(const std::string &filepath,
                 const std::vector<BinaryVector> &points) {
  if (points.empty()) {
    throw std::runtime_error("Cannot write empty points set.");
  }

  Index dimPoints = points[0].size();
  for (const auto &point : points) {
    if (point.size() != dimPoints) {
      throw std::runtime_error("Inconsistent point dimensions.");
    }
  }

  std::ofstream outfile(filepath);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open file for writing.");
  }

  // Write header
  outfile << "points " << dimPoints << " " << points.size() << "\n";

  // Write points
  for (const auto &point : points) {
    for (Index i = 0; i < point.size(); ++i) {
      outfile << static_cast<int>(point[i]);
      if (i < point.size() - 1)
        outfile << " ";
    }
    outfile << "\n";
  }
}

// ============================================================================
// Domain I/O
// ============================================================================

Domain readDomain(const std::string &filepath) {
  std::ifstream infile(filepath);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file.");
  }

  std::string header;
  if (!std::getline(infile, header)) {
    throw std::runtime_error("File is empty.");
  }

  std::istringstream ss(header);
  std::string keyword;
  Index dimConstr;
  Index numConstr;

  if (!(ss >> keyword >> dimConstr >> numConstr)) {
    throw std::runtime_error("Invalid header format.");
  } else if (keyword != "domain") {
    throw std::runtime_error("Invalid header format.");
  }

  if (dimConstr <= 0) {
    throw std::runtime_error("Constraint dimension must be positive.");
  }

  Domain domain;
  domain.reserve(numConstr);

  std::string line;
  Index constrCount = 0;
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;

    std::istringstream row_stream(line);
    RealVector coeffs;
    coeffs.reserve(dimConstr);

    // Read the dimConstr coefficients for vector 'a'
    double val;
    for (Index i = 0; i < dimConstr; ++i) {
      if (!(row_stream >> val)) {
        throw std::runtime_error("Not enough coefficients in constraint row.");
      }
      coeffs.push_back(val);
    }

    // Read the scalar 'b'
    double b;
    if (!(row_stream >> b)) {
      throw std::runtime_error("Missing scalar b in constraint row.");
    }

    // Read the relation symbol
    std::string relSymbol;
    if (!(row_stream >> relSymbol)) {
      throw std::runtime_error("Missing relation symbol in constraint row.");
    }

    Relation relation = stringToRelationType(relSymbol);

    // Fill the constraint tuple and add to domain
    domain.emplace_back(std::move(coeffs), b, relation);
    constrCount++;
  }

  if (constrCount != numConstr) {
    throw std::runtime_error("Incorrect number of constraints.");
  }

  return domain;
}

void writeDomain(const std::string &filepath, const Domain &domain) {
  if (domain.empty()) {
    throw std::runtime_error("Cannot write empty domain.");
  }

  Index dimConstr = std::get<0>(domain[0]).size();
  for (const auto &constr : domain) {
    if (std::get<0>(constr).size() != dimConstr) {
      throw std::runtime_error("Inconsistent constraint dimensions.");
    }
  }

  std::ofstream outfile(filepath);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open file for writing.");
  }

  // Write header
  outfile << "domain " << dimConstr << " " << domain.size() << "\n";

  // Write constraints
  outfile << std::setprecision(17);
  for (const auto &constr : domain) {

    const auto &[a, b, r] = constr;

    // Write coefficients
    for (Index i = 0; i < a.size(); ++i) {
      outfile << a[i];
      outfile << " ";
    }

    // Write scalar b
    outfile << b << " ";

    // Write relation symbol
    outfile << relationTypeToString(r);

    outfile << "\n";
  }
}

// ============================================================================
// Queries I/O
// ============================================================================

std::vector<RealVector> readQueries(const std::string &filepath) {
  std::ifstream infile(filepath);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file.");
  }

  std::string header;
  if (!std::getline(infile, header)) {
    throw std::runtime_error("File is empty.");
  }

  std::istringstream ss(header);
  std::string keyword;
  Index dimQueries;
  Index numQueries;

  if (!(ss >> keyword >> dimQueries >> numQueries)) {
    throw std::runtime_error("Invalid header format.");
  } else if (keyword != "queries") {
    throw std::runtime_error("Invalid header format.");
  }

  if (numQueries <= 0 || dimQueries <= 0) {
    throw std::runtime_error(
        "Number and dimension of queries must be positive.");
  }

  std::vector<RealVector> queries;
  queries.reserve(numQueries);

  std::string line;
  Index queriesCount = 0;
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;

    std::istringstream row_stream(line);
    RealVector query;
    query.reserve(dimQueries);

    double val;
    while (row_stream >> val) {
      query.push_back(val);
    }

    if (query.size() != dimQueries) {
      throw std::runtime_error("Incorrect query dimension.");
    }

    queries.push_back(std::move(query));
    queriesCount++;
  }

  if (queriesCount != numQueries) {
    throw std::runtime_error("Incorrect number of queries.");
  }

  return queries;
}

void writeQueries(const std::string &filepath,
                  const std::vector<RealVector> &queries) {
  if (queries.empty()) {
    throw std::runtime_error("Cannot write empty queries set.");
  }

  Index dimQueries = queries[0].size();
  for (const auto &query : queries) {
    if (query.size() != dimQueries) {
      throw std::runtime_error("Inconsistent query dimensions.");
    }
  }

  std::ofstream outfile(filepath);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open file for writing.");
  }

  // Write header
  outfile << "queries " << dimQueries << " " << queries.size() << "\n";

  // Write queries
  outfile << std::setprecision(17);
  for (const auto &query : queries) {
    for (Index i = 0; i < query.size(); ++i) {
      outfile << query[i];
      if (i < query.size() - 1)
        outfile << " ";
    }
    outfile << "\n";
  }
}

} // namespace treeco