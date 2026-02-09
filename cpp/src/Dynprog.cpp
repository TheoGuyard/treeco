#include "treeco/Dynprog.hpp"

namespace treeco {

std::ostream& operator<<(std::ostream& oss, const DynprogStats& stats) {
  oss << "Dynprog statistics\n";
  oss << "  total time     : " << std::fixed << std::setprecision(4) << stats.runTime << "\n";
  oss << "  iterations     : " << stats.numIters << "\n";
  oss << "  evaluations    : " << stats.numEvals << "\n";
  oss << "  states created : " << stats.numStates << "\n";
  oss << "  states closed  : " << stats.numStatesClosed << "\n";
  oss << "  states pruned  : " << stats.numStatesPruned << "\n";
  oss << "  lp solved      : " << stats.lpSolved << "\n";
  if (stats.optimalDepth == MAX_DEPTH) {
    oss << "  optimal depth  : inf (not found)";
  } else {
    oss << "  optimal depth  : " << stats.optimalDepth << "";
  }
  return oss;
}

Dynprog::Dynprog(const Voronoi& voronoi, const Domain& domain) : voronoi_(voronoi), domain_(domain) {
  if (!voronoi_.isBuilt()) { throw std::invalid_argument("Voronoi diagram must be built before using Dynprog"); }
}

void Dynprog::run(const DynprogParams& params) {
  // ------------------------------------------
  // Setup attributes
  // ------------------------------------------

  params_ = params;
  stats_ = DynprogStats();
  status_ = DynprogStatus::INVALID;
  startTime_ = Clock::now();
  checkTime_ = startTime_;

  switch (params_.branching) {
    case Branching::BINARY:
      branchDirections_ = {Relation::LT, Relation::GT};
      break;
    case Branching::TERNARY:
      branchDirections_ = {Relation::LT, Relation::GT, Relation::EQ};
      break;
    default:
      throw std::invalid_argument("Unknown branching strategy");
  }

  // ------------------------------------------
  // Initialization
  // ------------------------------------------

  logHeader();

  // Initialize feasibility checker with splits and input domain specs
  feasibility_ = std::make_unique<Feasibility>(voronoi_.splits(), domain_, params_.tolerance, params_.useSlacks);

  states_.clear();
  regionToStateId_.clear();
  initPositions();
  initRootState();

  // ------------------------------------------
  // Iterative exploration loop
  // ------------------------------------------

  logProgress("initialization");
  logSave();

  for (Index kSplits : getIterationRange()) {
    stats_.numIters = kSplits;

    evaluateState(rootId_, kSplits);
    updateStatus();

    if (status_ == DynprogStatus::OPTIMAL) {
      logProgress("optimal");
      logSave();
      break;
    } else if (checkTimeLimit()) {
      logProgress("time limit");
      logSave();
      break;
    }

    logProgress("subopt (" + std::to_string(kSplits) + " splits)");
    logSave();
  }

  // ------------------------------------------
  // Finalization
  // ------------------------------------------

  // Fill remaining statistics fields
  stats_.runTime = std::chrono::duration<double>(Clock::now() - startTime_).count();
  stats_.numStates = states_.size();
  stats_.lpSolved = feasibility_->numSolve();
  stats_.optimalDepth = states_[rootId_].ubHeight;

  logFooter();

  return;
}

void Dynprog::initPositions() {
  // Initialize position cache container with RF values (unknown position)
  positions_.resize(voronoi_.numFaces());
  for (auto& facePositions : positions_) { facePositions.resize(voronoi_.numSplits(), Relation::RF); }

  // Precompute all positions
  for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) {
#ifdef TREECO_BUILD_PYTHON
    // Propagate keyboard interrupts from python
    if (Py_IsInitialized()) {
      py::gil_scoped_acquire acquire;
      if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
    }
#endif

    const Cone& face = voronoi_.face(faceId).cone;
    const Cone& faceInt = face.interior();

    feasibility_->add(face.cuts());
    for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) {
      // Check if the split if a face boundary
      for (const auto& [sid, dir] : face.cuts()) {
        if (sid == splitId) {
          positions_[faceId][splitId] = dir;
          break;
        }
      }

      if (positions_[faceId][splitId] != Relation::RF) { continue; }

      // Determine position with feasibility check if not a boundary
      // - face ∩ s> = ∅ → face ⊆ s<=
      //    - int(face) ∩ s = ∅ → face ⊆ s< ∪ {0}
      //    - int(face) ∩ s ≠ ∅ → face ⊆ s< ∪ {s=}
      // - face ∩ s< = ∅ → face ⊆ s>=
      //    - int(face) ∩ s = ∅ → face ⊆ s> ∪ {0}
      //    - int(face) ∩ s ≠ ∅ → face ⊆ s> ∪ {s=}
      // - otherwise → int(face) intersects s< and s>

      Cut cutGT = Cut(splitId, Relation::GT);
      feasibility_->add(cutGT);
      bool checkGT = feasibility_->check();
      feasibility_->remove(cutGT);
      if (!checkGT) {
        feasibility_->add(faceInt.cuts());
        bool checkGE = feasibility_->check();
        feasibility_->remove(faceInt.cuts());
        if (!checkGE) {
          positions_[faceId][splitId] = Relation::LT;
        } else {
          positions_[faceId][splitId] = Relation::LE;
        }
      } else {
        Cut cutLT = Cut(splitId, Relation::LT);
        feasibility_->add(cutLT);
        bool checkLT = feasibility_->check();
        feasibility_->remove(cutLT);
        if (!checkLT) {
          feasibility_->add(faceInt.cuts());
          bool checkLE = feasibility_->check();
          feasibility_->remove(faceInt.cuts());
          if (!checkLE) {
            positions_[faceId][splitId] = Relation::GT;
          } else {
            positions_[faceId][splitId] = Relation::GE;
          }
        } else {
          positions_[faceId][splitId] = Relation::EQ;
        }
      }
    }
    feasibility_->remove(face.cuts());
  }
}

void Dynprog::initRootState() {
  // Initialize root state
  State root = State{Cone()};

  // Initialize root faces (all faces if no input domain, else filtered)
  root.faceIds.reserve(voronoi_.numFaces());
  if (domain_.empty()) {
    for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) { root.faceIds.push_back(faceId); }
  } else {
    for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) {
      const Face& face = voronoi_.face(faceId);
      feasibility_->add(face.cone.cuts());
      if (feasibility_->check()) { root.faceIds.push_back(faceId); }
      feasibility_->remove(face.cone.cuts());
    }
  }
  root.faceIds.shrink_to_fit();

  // Initialize root splits (all unbuild splits)
  root.splits.reserve(voronoi_.numSplits());
  for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) { root.splits.push_back(Split{splitId}); }

  // Initialize root lower bound
  evaluateLowerBound(root);

  // Register root state
  rootId_ = states_.size();
  regionToStateId_[root.region] = rootId_;
  states_.push_back(root);
}

Index Dynprog::getOrCreateChildId(Index parentId, Index splitPos, Relation cutDir) {
  State& parent = states_[parentId];
  Split& split = states_[parentId].splits[splitPos];

  const Cut childCut = Cut(split.id, cutDir);
  const Cone childCone = parent.region.refine(childCut);

  // Recover child id from memoization or create new child
  Index childId = INVALID_INDEX;
  auto itChild = regionToStateId_.find(childCone);
  if (itChild != regionToStateId_.end()) {
    childId = itChild->second;
  } else {
    State child = createState(parent, split, cutDir);
    childId = states_.size();
    regionToStateId_[child.region] = childId;
    states_.push_back(std::move(child));
  }

  return childId;
}

State Dynprog::createState(const State& parent, const Split& split, Relation cutDir) {
  Cut cut = Cut(split.id, cutDir);
  Cone region = parent.region.refine(cut);
  State child = State{region};

  // Initialize child faces
  child.faceIds.reserve(parent.faceIds.size());
  feasibility_->add(child.region.cuts());
  for (Index faceId : parent.faceIds) {
    bool intersect = checkChildFace(parent, split, cutDir, faceId, true);
    if (intersect) { child.faceIds.push_back(faceId); }
  }
  feasibility_->remove(child.region.cuts());
  child.faceIds.shrink_to_fit();

  // Initialize child splits
  child.splits.reserve(std::max(parent.splits.size() - 1, Index(0)));
  for (const Split& parentSplit : parent.splits) {
    if (parentSplit.id == split.id) { continue; }
    child.splits.push_back(Split{parentSplit.id});
  }
  child.splits.shrink_to_fit();

  // Initialize child lower bound
  evaluateLowerBound(child);

  // Test if new child is a leaf
  if (child.isLeaf()) {
    child.lbHeight = 0;
    child.ubHeight = 0;
    child.isClosed = true;
    stats_.numStatesLeafed++;
  }

  return child;
}

State Dynprog::dummyChild(const State& parent, const Split& split, Relation cutDir, const Cone& region) {
  State child = State{region};

  child.faceIds.reserve(parent.faceIds.size());
  for (Index faceId : parent.faceIds) {
    std::optional<bool> isEmpty = inferChildFace(parent, split, cutDir, faceId);
    if (isEmpty.has_value() && isEmpty.value()) { child.faceIds.push_back(faceId); }
  }
  child.faceIds.shrink_to_fit();

  evaluateLowerBound(child);

  return child;
}

void Dynprog::buildState(Index stateId) {
  stats_.numStatesBuilt++;

  State& state = states_[stateId];
  std::vector<Split>& splits = state.splits;

  // ----- Build splits sequentially ----- //

  // Buffer to store new children created newChildren[splitPos][cutDir] = child
  std::map<Index, std::map<Index, State>> newChildren = {};

  double bestSplitScore = std::numeric_limits<double>::infinity();
  Index bestSplitId = INVALID_INDEX;
  Index bestSplitLb = MAX_DEPTH;
  std::vector<Index> invalidSplitIds = {};

  feasibility_->add(state.region.cuts());
  for (Index splitPos = 0; splitPos < splits.size(); ++splitPos) {
    Split& split = splits[splitPos];

    // Split child data
    Index childNonempty = 0;
    Index childIdentical = 0;
    Index childBestLb = 0;
    std::vector<Index> childFaceCounts = {};

    for (Index dirId = 0; dirId < branchDirections_.size(); ++dirId) {
      Relation cutDir = branchDirections_[dirId];
      Cut childCut = Cut(split.id, cutDir);
      Cone childRegion = state.region.refine(childCut);

      // Recover child indices from memoization
      // If not found and strong checks enabled: instantiate new child state
      auto itChild = regionToStateId_.find(childRegion);
      if (itChild != regionToStateId_.end()) {
        split.childIds[dirId] = itChild->second;
      } else {
        split.childIds[dirId] = INVALID_INDEX;
        if (params_.strongChecks) {
          State child = createState(state, split, cutDir);
          newChildren[splitPos][dirId] = std::move(child);
        }
      }

      // Recover child (dummy child if not found and strong checks disabled)
      Index childId = split.childIds[dirId];
      const State& child = (childId != INVALID_INDEX) ? states_[childId]
                           : params_.strongChecks     ? newChildren[splitPos][dirId]
                                                      : dummyChild(state, split, cutDir, childRegion);

      // Update child data
      if (child.faceIds.size() > 0) { childNonempty += 1; }
      if (child.faceIds == state.faceIds) { childIdentical += 1; }
      childBestLb = std::max(childBestLb, safeAdd(child.lbHeight));
      childFaceCounts.push_back(child.faceIds.size());

      // Infer split validity
      if (params_.filterChecks) {
        if (childNonempty >= 2) { split.valid = true; }
        if (childIdentical >= 1) { split.valid = false; }
      }

      // Break if split found invalid
      if (split.valid.has_value() && !split.valid.value()) { break; }
    }

    // Check split validity if not already inferred
    if (!split.valid.has_value()) {
      Cut cut = Cut(split.id, Relation::EQ);
      feasibility_->add(state.region.cuts());
      feasibility_->add(cut);
      split.valid = feasibility_->check();
      feasibility_->remove(cut);
      feasibility_->remove(state.region.cuts());
    }

    // Mark invalid splits and continue
    if (!split.valid.value()) {
      invalidSplitIds.push_back(split.id);
      continue;
    }

    // Evaluate split score
    evaluateScore(state, split, childFaceCounts);

    // Update best split score, id, and lower bound
    if (split.score < bestSplitScore) {
      bestSplitScore = split.score;
      bestSplitId = split.id;
    }
    if (childBestLb < bestSplitLb) { bestSplitLb = childBestLb; }
  }
  feasibility_->remove(state.region.cuts());

  // ----- Insert new children in memoization ----- //

  // Warning: state and splits references are dangling from here

  for (auto& [splitPos, splitChildren] : newChildren) {
    for (auto& [dirId, child] : splitChildren) {
      // For GREEDY exploration, only instantiate children of the best split
      if (params_.exploration == Exploration::GREEDY) {
        if (states_[stateId].splits[splitPos].id != bestSplitId) { continue; }
      }

      // Skip children of invalid splits
      if (!states_[stateId].splits[splitPos].valid.value()) { continue; }

      // Remove any pending invalid child split
      child.splits.erase(std::remove_if(child.splits.begin(), child.splits.end(),
                                        [&](const Split& childSplit) {
                                          auto it = find(invalidSplitIds.begin(), invalidSplitIds.end(), childSplit.id);
                                          return it != invalidSplitIds.end();
                                        }),
                         child.splits.end());

      // Register child in memoization
      Index childId = states_.size();
      regionToStateId_[child.region] = childId;
      states_[stateId].splits[splitPos].childIds[dirId] = childId;
      states_.push_back(std::move(child));
    }
  }

  // ----- Register building results in state ----- //

  // Register that the state has been built
  states_[stateId].isBuilt = true;

  // Remove invalid splits
  states_[stateId].splits.erase(std::remove_if(std::begin(states_[stateId].splits), std::end(states_[stateId].splits),
                                               [](const auto& split) { return !split.valid.value(); }),
                                std::end(states_[stateId].splits));

  // Check if the state became a leaf (no more valid splits)
  if (states_[stateId].splits.size() == 0) {
    states_[stateId].lbHeight = 0;
    states_[stateId].ubHeight = 0;
    states_[stateId].isClosed = true;
    stats_.numStatesLeafed++;
    return;
  }

  // Sort splits by ascending score
  std::stable_sort(std::begin(states_[stateId].splits), std::end(states_[stateId].splits),
                   [](const Split& splitA, const Split& splitB) { return splitA.score < splitB.score; });

  // Update the state lower bound
  if (bestSplitLb != MAX_DEPTH) { states_[stateId].lbHeight = std::max(states_[stateId].lbHeight, bestSplitLb); }

  // Check if the state can be pruned with updated lower bound
  if (states_[stateId].depth() + states_[stateId].lbHeight >= states_[rootId_].ubHeight) {
    states_[stateId].lbHeight = MAX_DEPTH;
    states_[stateId].ubHeight = MAX_DEPTH;
    states_[stateId].isClosed = true;
    stats_.numStatesPruned++;
  }
}

void Dynprog::evaluateState(Index stateId, Index kSplits) {
#ifdef TREECO_BUILD_PYTHON
  // Propagate keyboard interrupts from python
  if (Py_IsInitialized()) {
    py::gil_scoped_acquire acquire;
    if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
  }
#endif

  // Break exploration if time limit exceeded
  if (checkTimeLimit()) { return; }

  // -------------------- State pre-processing -------------------- //

  stats_.numEvals++;

  // Return if state is closed (leaf, pruned, or optimal)
  if (states_[stateId].isClosed) { return; }

  // Check if the state is a leaf
  if (states_[stateId].isLeaf()) {
    states_[stateId].lbHeight = 0;
    states_[stateId].ubHeight = 0;
    states_[stateId].isClosed = true;
    stats_.numStatesLeafed++;
    return;
  }

  // Check if the state can be pruned
  if (states_[stateId].depth() + states_[stateId].lbHeight >= states_[rootId_].ubHeight) {
    states_[stateId].lbHeight = MAX_DEPTH;
    states_[stateId].ubHeight = MAX_DEPTH;
    states_[stateId].isClosed = true;
    stats_.numStatesPruned++;
    return;
  }

  // Check if the state splits need to be built (and if closed while building)
  if (!states_[stateId].isBuilt) {
    buildState(stateId);
    if (states_[stateId].isClosed) { return; }
  }

  // -------------------- State processing -------------------- //

  Index bestSplitLb = MAX_DEPTH;
  for (Index splitPos = 0; splitPos < states_[stateId].splits.size(); ++splitPos) {
    // Break if enough splits have been tested
    if (splitPos >= kSplits) { break; }

    logProgress();

    // Skip if the split is closed (i.e., pruned, or optimal)
    if (states_[stateId].splits[splitPos].isClosed) { continue; }

    // Process each children
    Index lbHeightSplit = 0;
    Index ubHeightSplit = 0;
    for (Index i = 0; i < branchDirections_.size(); ++i) {
      Index childId = states_[stateId].splits[splitPos].childIds[i];

      // Recover child id if needed
      if (childId == INVALID_INDEX) {
        childId = getOrCreateChildId(stateId, splitPos, branchDirections_[i]);
        states_[stateId].splits[splitPos].childIds[i] = childId;
      }

      // Update the split lower bound
      lbHeightSplit = std::max(lbHeightSplit, safeAdd(states_[childId].lbHeight));

      // Check if the split can be pruned
      if (states_[stateId].depth() + lbHeightSplit >= states_[0].ubHeight) {
        states_[stateId].splits[splitPos].isClosed = true;
        break;
      }

      // Process the child
      evaluateState(childId, kSplits);

      // Update the split upper bound
      ubHeightSplit = std::max(ubHeightSplit, safeAdd(states_[childId].ubHeight));

      // Break if the split cannot provide any improvement
      if (ubHeightSplit >= states_[stateId].ubHeight) { break; }
    }

    // Skip next steps if the split was closed during children processing
    if (states_[stateId].splits[splitPos].isClosed) { continue; }

    // Update best split lower bound
    bestSplitLb = std::min(bestSplitLb, lbHeightSplit);

    // Update the best state upper bound and split index
    if (ubHeightSplit < states_[stateId].ubHeight) {
      states_[stateId].splitId = states_[stateId].splits[splitPos].id;
      states_[stateId].ubHeight = ubHeightSplit;
    }

    // Check if the split is optimal by matching bounds
    if (lbHeightSplit == ubHeightSplit) { states_[stateId].splits[splitPos].isClosed = true; }
  }

  // -------------------- State post-processing -------------------- //

  // Check if the state is optimal by matching bounds
  if (states_[stateId].lbHeight == states_[stateId].ubHeight) {
    states_[stateId].isClosed = true;
    stats_.numStatesClosed++;
    return;
  }

  // Check if the state is optimal by splits exhaustion
  if (kSplits >= states_[stateId].splits.size()) {
    states_[stateId].lbHeight = states_[stateId].ubHeight;
    states_[stateId].isClosed = true;
    stats_.numStatesClosed++;
  }
}

void Dynprog::evaluateLowerBound(State& state) {
  Index numFaces = state.faceIds.size();

  if (numFaces <= 1) {
    state.lbHeight = 0;
    return;
  }

  Index base = branchDirections_.size();

  state.lbHeight =
      static_cast<Index>(std::ceil(std::log(static_cast<double>(numFaces)) / std::log(static_cast<double>(base))));
}

void Dynprog::evaluateScore(State& state, Split& split, const std::vector<Index>& childFaceCounts) {
  Index stateFaceCount = state.faceIds.size();
  Index numChildren = childFaceCounts.size();

  double score;

  switch (params_.splitScoring) {
    // Variance from equal split: lower is better
    case SplitScoring::VARIANCE: {
      double mean = static_cast<double>(stateFaceCount) / numChildren;
      double svar = 0.0;
      for (Index count : childFaceCounts) {
        double diff = static_cast<double>(count) - mean;
        svar += diff * diff;
      }
      score = svar / numChildren;
      break;
    }

    // Information gain: higher reduction is better (set negative)
    case SplitScoring::ENTROPY: {
      Index totalChildFaces = 0;
      for (Index count : childFaceCounts) { totalChildFaces += count; }
      if (totalChildFaces == 0) {
        score = std::numeric_limits<double>::infinity();
        break;
      }
      double stateEntropy = std::log2(static_cast<double>(stateFaceCount));
      double childEntropy = 0.0;
      for (Index count : childFaceCounts) {
        if (count > 0) {
          double p = static_cast<double>(count) / totalChildFaces;
          childEntropy += p * std::log2(static_cast<double>(count));
        }
      }
      double gainEntropy = stateEntropy - childEntropy;
      score = -gainEntropy;
      break;
    }

    // Minimize maximum child height: lower is better
    case SplitScoring::MINMAX: {
      Index maxCount = 0;
      for (Index count : childFaceCounts) { maxCount = std::max(maxCount, count); }
      score = static_cast<double>(maxCount);
      break;
    }

    // Random scoring in [0,1)
    case SplitScoring::RANDOM: {
      std::mt19937 rng_;
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      score = dist(rng_);
      break;
    }

    // No scoring: all splits equal
    case SplitScoring::NONE: {
      score = 0.0;
      break;
    }

    default:
      throw std::invalid_argument("Unknown split scoring strategy");
  }

  split.score = score;
}

bool Dynprog::checkChildFace(const State& state, const Split& split, Relation cutDir, Index faceId,
                             bool externalChildCuts) {
  // Try to infer whether the face is a child one based on filtering rules
  if (params_.filterChecks) {
    std::optional<bool> isInferredChildFace = inferChildFace(state, split, cutDir, faceId);
    if (isInferredChildFace.has_value()) { return isInferredChildFace.value(); }
  }

  // If filtering is not enabled or failed, perform full feasibility check
  if (!externalChildCuts) {
    feasibility_->add(state.region.cuts());
    feasibility_->add(Cut(split.id, cutDir));
  }
  feasibility_->add(voronoi_.face(faceId).cone.cuts());
  bool isChildFace = feasibility_->check();
  feasibility_->remove(voronoi_.face(faceId).cone.cuts());
  if (!externalChildCuts) {
    feasibility_->remove(state.region.cuts());
    feasibility_->remove(Cut(split.id, cutDir));
  }

  return isChildFace;
}

void Dynprog::updateStatus() {
  if (states_[rootId_].isClosed) {
    status_ = DynprogStatus::OPTIMAL;
  } else if (states_[rootId_].ubHeight < MAX_DEPTH) {
    status_ = DynprogStatus::SUBOPTIMAL;
  } else {
    status_ = DynprogStatus::INVALID;
  }
}

std::optional<bool> Dynprog::inferChildFace(const State& state, const Split& split, Relation cutDir, Index faceId) {
  Relation position = positions_[faceId][split.id];
  const Cone& stateCone = state.region;

  if (cutDir == Relation::LT) {
    if (position == Relation::GE || position == Relation::GT) {
      return false;
    } else if (position == Relation::LE || position == Relation::LT) {
      if (stateCone.isOpen()) { return true; }
    } else if (position == Relation::EQ) {
      // no inference possible
    } else {
      throw std::invalid_argument("Unknown relation type");
    }
  } else if (cutDir == Relation::GT) {
    if (position == Relation::LE || position == Relation::LT) {
      return false;
    } else if (position == Relation::GE || position == Relation::GT) {
      if (stateCone.isOpen()) { return true; }
    } else if (position == Relation::EQ) {
      // no inference possible
    } else {
      throw std::invalid_argument("Unknown relation type");
    }
  } else if (cutDir == Relation::EQ) {
    if (stateCone.isOpen()) {
      if (position == Relation::LT || position == Relation::GT) { return false; }
    }
  } else {
    throw std::invalid_argument("Unknown cut direction");
  }

  // Ultimate check if no inference passed
  if (childContainsFaceCenter(state, split, cutDir, faceId)) { return true; }

  // Fallback: no inference possible
  return std::nullopt;
}

bool Dynprog::childContainsFaceCenter(const State& state, const Split& split, Relation cutDir, Index faceId) {
  const Cone& childCone = state.region.refine(Cut(split.id, cutDir));
  const SimplexVector& faceCenter = voronoi_.point(faceId);

  for (const auto& [splitId, cutDir] : childCone.cuts()) {
    int scalProd = dot(voronoi_.split(splitId), faceCenter);
    switch (cutDir) {
      case Relation::LT:
        if (scalProd >= -params_.tolerance) { return false; }
        break;
      case Relation::GT:
        if (scalProd <= params_.tolerance) { return false; }
        break;
      case Relation::EQ:
        if (std::abs(scalProd) > params_.tolerance) { return false; }
        break;
      default:
        throw std::invalid_argument("Unknown relation type");
    }
  }

  return true;
}

std::vector<Index> Dynprog::getIterationRange() const {
  std::vector<Index> iterationRange;
  switch (params_.exploration) {
    case Exploration::GREEDY:
      iterationRange.push_back(1);
      break;
    case Exploration::ITERATIVE:
      iterationRange.reserve(voronoi_.numSplits());
      for (Index k = 1; k <= voronoi_.numSplits(); ++k) { iterationRange.push_back(k); }
      break;
    case Exploration::EXHAUSTIVE:
      iterationRange.push_back(voronoi_.numSplits());
      break;
    default:
      throw std::invalid_argument("Unknown exploration strategy");
  }
  return iterationRange;
}

bool Dynprog::checkTimeLimit() const {
  double elapsed = std::chrono::duration<double>(Clock::now() - startTime_).count();
  return elapsed >= params_.timeLimit;
}

void Dynprog::logHeader() const {
  if (!params_.verbose) return;
  std::ostream& out = *(params_.outputStream);
  out << "Running dynamic programming...\n";
  out << "  num faces  : " << voronoi_.numFaces() << "\n";
  out << "  num splits : " << voronoi_.numSplits() << "\n";
  out << std::string("  ") + std::string(9 * 12, '-') << "\n";
  out << "  ";
  out << std::setw(12) << "time";
  out << std::setw(12) << "iters";
  out << std::setw(12) << "lb depth";
  out << std::setw(12) << "ub depth";
  out << std::setw(12) << "lp opt";
  out << std::setw(12) << "states";
  out << std::setw(12) << "built";
  out << std::setw(12) << "closed";
  out << std::setw(12) << "evals";
  out << "\n";
}

void Dynprog::logProgress(const std::string& message) {
  if (!params_.verbose) { return; }
  if (elapsedTime(checkTime_) < params_.logInterval && message.empty()) { return; }

  checkTime_ = Clock::now();

  std::ostream& out = *(params_.outputStream);
  out << "  ";
  out << std::setw(12) << std::fixed << std::setprecision(2) << elapsedTime(startTime_);
  out << std::setw(12) << stats_.numIters;
  out << std::setw(12) << states_[rootId_].lbHeight;
  out << std::setw(12) << (states_[rootId_].ubHeight == MAX_DEPTH ? "inf" : std::to_string(states_[rootId_].ubHeight));
  out << std::setw(12) << feasibility_->numSolve();
  out << std::setw(12) << states_.size();
  out << std::setw(12) << stats_.numStatesClosed;
  out << std::setw(12) << stats_.numStatesBuilt;
  out << std::setw(12) << stats_.numEvals;
  out << "  " << message;
  out << "\n";
}

void Dynprog::logSave() {
  if (!params_.logSave) return;
  stats_.runTime = elapsedTime(startTime_);
  stats_.numStates = states_.size();
  stats_.lpSolved = feasibility_->numSolve();
  stats_.optimalDepth = states_[rootId_].ubHeight;
  logs_.push_back(stats_);
}

void Dynprog::logFooter() const {
  if (!params_.verbose) return;
  std::ostream& out = *(params_.outputStream);
  out << std::string("  ") + std::string(9 * 12, '-') << "\n";
  out << "  time: " << std::fixed << std::setprecision(4) << stats_.runTime << "\n";
}

}  // namespace treeco
