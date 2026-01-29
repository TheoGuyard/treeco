#include "treeco/Tree.hpp"

namespace treeco {

Tree::Tree(const Voronoi& voronoi) : voronoi_(voronoi) {}

void Tree::synthetize(const Dynprog& dynprog, const TreeParams& params) {
  clear();

  params_ = params;
  startTime_ = Clock::now();
  checkTime_ = startTime_;

  logHeader();

  if (dynprog.status() == DynprogStatus::INVALID) {
    throw std::runtime_error("Dynprog terminate with an invalid status");
  }

  std::queue<std::pair<Index, Index>> queue;
  Index stateId = dynprog.rootId();
  const State& state = dynprog.state(stateId);
  Index nodeId =
      addRoot(state.isLeaf() ? state.faceIds : std::vector<Index>{}, state.isLeaf() ? INVALID_INDEX : state.splitId);
  queue.emplace(stateId, nodeId);

  for (; !queue.empty(); queue.pop()) {
    auto [stateId, nodeId] = queue.front();
    const Node& node = nodes_.at(nodeId);
    if (node.type == NodeType::NODE) {
      const State& state = dynprog.state(stateId);
      for (const auto& [relation, childStateId] : state.splits.at(state.splitId).childIds) {
        const State& childState = dynprog.state(childStateId);
        Index childNodeId = addNode(nodeId, relation, childState.isLeaf() ? childState.faceIds : std::vector<Index>{},
                                    childState.isLeaf() ? INVALID_INDEX : childState.splitId);
        queue.emplace(childStateId, childNodeId);
      }
    }
    logProgress();

#ifdef TREECO_BUILD_PYTHON
    // Propagate keyboard interrupts from python
    if (Py_IsInitialized()) {
      py::gil_scoped_acquire acquire;
      if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
    }
#endif
  }

  logProgress("done");

  stats_.isBuilt = true;
  stats_.buildTime = elapsedTime(startTime_);
  stats_.dynprogStats = dynprog.stats();
  stats_.dynprogLogs = dynprog.logs();

  logFooter();
}

Index Tree::addRoot(const std::vector<Index>& pointsIds, Index splitId) {
  if (nodes_.size() > 0) { throw std::runtime_error("Tree is not empty"); }

  if (rootId_ != INVALID_INDEX) { throw std::runtime_error("Root already set"); }

  if (pointsIds.size() > 0 && splitId != INVALID_INDEX) {
    throw std::runtime_error("Cannot add root with both points and split");
  }

  NodeType rootType = splitId == INVALID_INDEX ? NodeType::LEAF : NodeType::NODE;

  Node root = Node{0, rootType, pointsIds, splitId};
  rootId_ = nodes_.size();
  nodes_.push_back(root);

  // Update tree attributes
  depth_ = 0;
  width_ = (root.type == NodeType::LEAF) ? 1 : 0;
  size_++;

  return rootId_;
}

Index Tree::addNode(Index parentId, Relation relation, const std::vector<Index>& pointsIds, Index splitId) {
  if (parentId < 0 || parentId >= nodes_.size()) { throw std::out_of_range("Parent index out of range"); }

  if (pointsIds.size() > 0 && splitId != INVALID_INDEX) {
    throw std::runtime_error("Cannot add node with both points and split");
  }

  Node& parent = nodes_.at(parentId);
  NodeType nodeType = splitId == INVALID_INDEX ? NodeType::LEAF : NodeType::NODE;

  if (parent.childIds.find(relation) != parent.childIds.end()) { throw std::runtime_error("Child already exists."); }

  Index nodeId = nodes_.size();

  Node node = Node{parent.depth + 1, nodeType, pointsIds, splitId};
  parent.childIds[relation] = nodeId;
  nodes_.push_back(node);

  // Update tree attributes
  depth_ = std::max(depth_, node.depth);
  width_ += (node.type == NodeType::LEAF) ? 1 : 0;
  size_++;

  return nodeId;
}

std::vector<SimplexVector> Tree::query(const std::vector<double>& cost) const {
  if (!isBuilt()) { throw std::runtime_error("Tree is not built yet."); }

  Index nodeId = rootId_;
  while (nodes_[nodeId].type == NodeType::NODE) {
    const Node& node = nodes_[nodeId];
    const TernaryVector& split = voronoi_.split(node.splitId);
    double splitSide = dot(split, cost);

    Relation relation;
    if (splitSide <= -params_.tolerance) {
      relation = Relation::LT;
    } else if (splitSide >= params_.tolerance) {
      relation = Relation::GT;
    } else {
      if (node.childIds.find(Relation::EQ) != node.childIds.end()) {
        relation = Relation::EQ;
      } else {
        relation = Relation::LT;
      }
    }
    nodeId = node.childIds.at(relation);
  }

  std::vector<SimplexVector> points;
  if (nodes_[nodeId].type == NodeType::LEAF) {
    points.reserve(nodes_[nodeId].pointsIds.size());
    for (int idx : nodes_[nodeId].pointsIds) { points.push_back(voronoi_.point(idx)); }
  } else {
    throw std::runtime_error("Reached a non-leaf node during query.");
  }

  return points;
}

void Tree::pprint(bool tightDisplay, std::ostream* outputStream) const {
  if (!isBuilt()) { throw std::runtime_error("Tree is not built yet."); }

  if (rootId_ == INVALID_INDEX) { throw std::runtime_error("Tree root is not defined."); }

  if (nodes_.empty()) {
    *params_.outputStream << "Empty tree\n";
    return;
  }

  pprintNode(nodes_.at(rootId_), "", "", true, tightDisplay, outputStream);
}

void Tree::pprintNode(const Node& node, const std::string& prefix, const std::string& label, bool last,
                      bool tightDisplay, std::ostream* outputStream) const {
  std::ostream& out = *outputStream;

  std::string branch = (label.empty() ? "" : (last ? "└── " : "├── "));
  out << prefix << branch << label;
  if (!label.empty()) { out << " "; }
  if (node.type == NodeType::LEAF) {
    out << "Leaf {";
    for (Index i = 0; i < node.pointsIds.size(); ++i) {
      if (tightDisplay) {
        out << node.pointsIds[i];
      } else {
        SimplexVector scaledPoint = voronoi_.point(node.pointsIds[i]);
        BinaryVector point = unscaleBinary(scaledPoint);
        printVector(point, outputStream);
      }
      if (i + 1 < node.pointsIds.size()) out << ",";
    }
    out << "}" << std::endl;
  } else {
    out << "Node (";
    if (tightDisplay) {
      out << node.splitId;
    } else {
      TernaryVector split = voronoi_.split(node.splitId);
      printVector(split, outputStream);
    }
    out << ")" << std::endl;
    std::string childPrefix = prefix + (label.empty() ? "" : (last ? "    " : "│   "));
    Index count = 0;
    for (const auto& [relation, childId] : node.childIds) {
      const Node& child = nodes_.at(childId);
      bool last = (++count == node.childIds.size());
      pprintNode(child, childPrefix, relationTypeToString(relation), last, tightDisplay, outputStream);
    }
  }
}

void Tree::clear() {
  nodes_.clear();
  stats_ = TreeStats();
  rootId_ = INVALID_INDEX;
  size_ = 0;
  width_ = 0;
  depth_ = 0;
}

void Tree::logHeader() const {
  if (!params_.verbose) return;
  std::ostream& out = *(params_.outputStream);
  out << "Synthesizing tree...\n";
  out << std::string("  ") + std::string(4 * 12, '-') << "\n";
  out << "  ";
  out << std::setw(12) << "time";
  out << std::setw(12) << "size";
  out << std::setw(12) << "width";
  out << std::setw(12) << "depth";
  out << "\n";
}

void Tree::logProgress(const std::string& message) {
  if (!params_.verbose) { return; }
  if (elapsedTime(checkTime_) < params_.logInterval && message.empty()) { return; }

  checkTime_ = Clock::now();

  std::ostream& out = *(params_.outputStream);
  out << "  ";
  out << std::setw(12) << std::fixed << std::setprecision(2) << elapsedTime(startTime_);
  out << std::setw(12) << size_;
  out << std::setw(12) << width_;
  out << std::setw(12) << depth_;
  out << "  " << message;
  out << "\n";
}

void Tree::logFooter() const {
  if (!params_.verbose) return;
  std::ostream& out = *(params_.outputStream);
  out << std::string("  ") + std::string(4 * 12, '-') << "\n";
  out << "  time: " << std::fixed << std::setprecision(4) << stats_.buildTime << "\n";
}

std::ostream& operator<<(std::ostream& oss, const Tree& tree) {
  oss << "Tree\n";
  oss << "  depth : " + std::to_string(tree.depth()) + "\n";
  oss << "  width : " + std::to_string(tree.width()) + "\n";
  oss << "  size  : " + std::to_string(tree.size()) + "\n";
  return oss;
}

}  // namespace treeco