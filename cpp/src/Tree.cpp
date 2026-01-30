#include "treeco/Tree.hpp"

namespace treeco {

void Tree::synthetize(const Dynprog& dynprog, const TreeParams& params) {
  clear();

  params_ = params;
  startTime_ = Clock::now();
  checkTime_ = startTime_;

  branchDirs_ = dynprog.branchDirections();

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

      auto splitIt = std::find_if(std::begin(state.splits), std::end(state.splits),
                                  [splitId = state.splitId](const Split& split) { return split.id == splitId; });

      for (Index i = 0; i < branchDirs_.size(); ++i) {
        Relation childDir = branchDirs_[i];
        Index childStateId = splitIt->childIds[i];

        const State& childState = dynprog.state(childStateId);
        Index childNodeId = addNode(nodeId, childDir, childState.isLeaf() ? childState.faceIds : std::vector<Index>{},
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

Index Tree::addNode(Index parentId, Relation childDir, const std::vector<Index>& pointsIds, Index splitId) {
  if (parentId < 0 || parentId >= nodes_.size()) { throw std::out_of_range("Parent index out of range"); }

  if (pointsIds.size() > 0 && splitId != INVALID_INDEX) {
    throw std::runtime_error("Cannot add node with both points and split");
  }

  Node& parent = nodes_.at(parentId);
  NodeType nodeType = splitId == INVALID_INDEX ? NodeType::LEAF : NodeType::NODE;

  for (const auto& child : parent.children) {
    if (child.first == childDir) { throw std::runtime_error("Parent already has a child with the same direction"); }
  }

  Index nodeId = nodes_.size();

  Node node = Node{parent.depth + 1, nodeType, pointsIds, splitId};
  parent.children.push_back(std::make_pair(childDir, nodeId));
  nodes_.push_back(node);

  // Update tree attributes
  depth_ = std::max(depth_, node.depth);
  width_ += (node.type == NodeType::LEAF) ? 1 : 0;
  size_++;

  return nodeId;
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