#include "treeco/Geometry.hpp"

namespace treeco {

// ----------------------------------------------------------------------------
// Cut implementation
// ----------------------------------------------------------------------------

Cut::Cut(Index hid, Relation dir) : hid(hid), dir(dir) {}

size_t Cut::hash() const {
  std::size_t h = std::hash<Index>{}(hid);
  h ^= std::hash<int>{}(static_cast<int>(dir)) + 0x9e3779b9 + (h << 6) + (h >> 2);
  return h;
}

std::string Cut::toString() const {
  std::ostringstream oss;
  oss << "Cut(" << hid << "," << relationTypeToString(dir) << ")";
  return oss.str();
}

bool Cut::operator==(const Cut& other) const { return hid == other.hid && dir == other.dir; }

bool Cut::operator!=(const Cut& other) const { return hid != other.hid || dir != other.dir; }

bool Cut::operator<(const Cut& other) const { return (hid < other.hid) || (hid == other.hid && dir < other.dir); }

// ----------------------------------------------------------------------------
// Cone implementation
// ----------------------------------------------------------------------------

Cone::Cone() = default;

Cone::Cone(const std::set<Cut>& cuts) {
  for (const auto& cut : cuts) { addCut(cut); }
}

Cone::Cone(const std::vector<Cut>& cuts) {
  for (const auto& cut : cuts) { addCut(cut); }
}

void Cone::addCut(const Cut& cut) {
  if (cuts_.size() == 0) { isOpen_ = false; }

  isOpen_ = isOpen_ ||
            (cut.dir == Relation::LT || cut.dir == Relation::GT || cut.dir == Relation::RF || cut.dir == Relation::RT);

  isPointed_ = isPointed_ && !(cut.dir == Relation::LT || cut.dir == Relation::GT || cut.dir == Relation::RF);

  cuts_.insert(cut);
}

void Cone::addCut(Index hid, Relation dir) { addCut(Cut(hid, dir)); }

const std::set<Cut>& Cone::cuts() const { return cuts_; }

Index Cone::numCuts() const { return cuts_.size(); }

bool Cone::isFullSpace() const { return cuts_.empty(); }

bool Cone::isOpen() const { return isOpen_; }

bool Cone::isPointed() const { return isPointed_; }

void Cone::clear() { cuts_.clear(); }

Cone Cone::refine(const Cut& cut) const {
  Cone refined = Cone(this->cuts_);
  refined.addCut(cut);
  return refined;
}

Cone Cone::interior() const {
  Cone interiorCone;
  for (const auto& [hid, dir] : cuts_) {
    Relation interiorDir = dir;
    switch (dir) {
      case Relation::LT:
        break;
      case Relation::LE:
        interiorDir = Relation::LT;
        break;
      case Relation::GT:
        break;
      case Relation::GE:
        interiorDir = Relation::GT;
        break;
      case Relation::EQ:
        interiorDir = Relation::RF;
        break;
      case Relation::RT:
        break;
      case Relation::RF:
        break;
      default:
        throw std::runtime_error("Invalid relation type.");
    }
    interiorCone.addCut(hid, interiorDir);
  }
  return interiorCone;
}

std::size_t Cone::hash() const {
  // Combine hashes of all halfspaces
  std::size_t h = 0;
  for (const auto& cut : cuts_) { h ^= cut.hash() + 0x9e3779b9 + (h << 6) + (h >> 2); }
  return h;
}

bool Cone::operator==(const Cone& other) const { return cuts_ == other.cuts(); }

std::ostream& operator<<(std::ostream& oss, const Cone& cone) {
  oss << "Cone{";
  for (const auto& cut : cone.cuts()) {
    oss << cut.toString();
    if (cut != *cone.cuts().rbegin()) { oss << ","; }
  }
  oss << "}";
  return oss;
}

// ----------------------------------------------------------------------------
// Utility functions
// ----------------------------------------------------------------------------

void normalize(TernaryVector& v) {
  for (size_t i = 0; i < v.size(); ++i) {
    if (v[i] != 0) {
      if (v[i] < 0) {
        for (size_t j = 0; j < v.size(); ++j) { v[j] = -v[j]; }
      }
      return;
    }
  }
}

int dot(const TernaryVector& p, const SimplexVector& v) {
  int result = 0;
  for (Index i = 0; i < p.size(); ++i) { result += p[i] * v[i]; }
  return result;
}

double dot(const TernaryVector& p, const RealVector& v) {
  double result = 0;
  for (Index i = 0; i < p.size(); ++i) { result += static_cast<double>(p[i]) * v[i]; }
  return result;
}

double dot(const RealVector& p, const BinaryVector& v) {
  double result = 0.0;
  for (Index i = 0; i < p.size(); ++i) { result += p[i] * static_cast<double>(v[i]); }
  return result;
}

double dot(const RealVector& p, const RealVector& v) {
  double result = 0.0;
  for (Index i = 0; i < p.size(); ++i) { result += p[i] * v[i]; }
  return result;
}

TernaryVector bisector(const SimplexVector& p1, const SimplexVector& p2) {
  TernaryVector normal;
  normal.resize(p1.size());
  for (size_t i = 0; i < p1.size(); ++i) { normal[i] = (p2[i] - p1[i]) / 2; }
  normalize(normal);
  return normal;
}

SimplexVector scaleBinary(const BinaryVector& x) {
  SimplexVector result(x.size());
  for (Index i = 0; i < x.size(); ++i) { result[i] = 2 * x[i] - 1; }
  return result;
}

BinaryVector unscaleBinary(const SimplexVector& x) {
  BinaryVector result(x.size());
  for (Index i = 0; i < x.size(); ++i) { result[i] = (x[i] + 1) / 2; }
  return result;
}

std::vector<SimplexVector> scaleBinarySet(const std::vector<BinaryVector>& X) {
  std::vector<SimplexVector> result;
  result.reserve(X.size());
  for (const auto& x : X) { result.push_back(scaleBinary(x)); }
  return result;
}

std::vector<BinaryVector> unscaleBinarySet(const std::vector<SimplexVector>& X) {
  std::vector<BinaryVector> result;
  result.reserve(X.size());
  for (const auto& x : X) { result.push_back(unscaleBinary(x)); }
  return result;
}

}  // namespace treeco
