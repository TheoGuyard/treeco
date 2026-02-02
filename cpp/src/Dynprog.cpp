#include "treeco/Dynprog.hpp"

namespace treeco {

std::ostream& operator<<(std::ostream& oss, const DynprogStats& stats) {
  oss << "Dynprog statistics\n";
  oss << "  total time     : " << std::fixed << std::setprecision(4)
      << stats.runTime << "\n";
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

Dynprog::Dynprog(const Voronoi& voronoi, const Domain& domain)
  : voronoi_(voronoi), domain_(domain) {
  if (!voronoi_.isBuilt()) {
    throw std::invalid_argument(
        "Voronoi diagram must be built before using Dynprog");
  }
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
  feasibility_ = std::make_unique<Feasibility>(voronoi_.splits(), domain_,
                                               params_.tolerance);

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
  stats_.runTime =
      std::chrono::duration<double>(Clock::now() - startTime_).count();
  stats_.numStates = states_.size();
  stats_.lpSolved = feasibility_->numSolve();
  stats_.optimalDepth = states_[rootId_].ubHeight;

  logFooter();

  return;
}

void Dynprog::initPositions() {
  // Initialize position cache container with RF values (unknown position)
  positions_.resize(voronoi_.numFaces());
  for (auto& facePositions : positions_) {
    facePositions.resize(voronoi_.numSplits(), Relation::RF);
  }

  // Precompute all positions
  for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) {
#ifdef TREECO_BUILD_PYTHON
    // Propagate keyboard interrupts from python
    if (Py_IsInitialized()) {
      py::gil_scoped_acquire acquire;
      if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
    }
#endif

    const Cone& interiorFace = voronoi_.face(faceId).cone.interior();
    feasibility_->add(interiorFace.cuts());
    for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) {
      const Cone& interiorCone = voronoi_.face(faceId).cone.interior();

      // Check if the split if a face boundary
      for (const auto& [sid, dir] : interiorCone.cuts()) {
        if (sid == splitId) {
          positions_[faceId][splitId] = dir;
          break;
        }
      }

      if (positions_[faceId][splitId] != Relation::RF) { continue; }

      // Determine position with feasibility check if not boundary
      // - int(face) ∩ s< = ∅ → face ⊆ s>=
      // - int(face) ∩ s> = ∅ → face ⊆ s<=
      // - otherwise → int(face) intersects s< and s>
      Cut cutLT = Cut(splitId, Relation::LT);
      feasibility_->add(cutLT);
      bool checkLT = feasibility_->check();
      feasibility_->remove(cutLT);
      if (!checkLT) {
        positions_[faceId][splitId] = Relation::GT;
      } else {
        Cut cutGT = Cut(splitId, Relation::GT);
        feasibility_->add(cutGT);
        bool checkGT = feasibility_->check();
        feasibility_->remove(cutGT);
        if (!checkGT) {
          positions_[faceId][splitId] = Relation::LT;
        } else {
          positions_[faceId][splitId] = Relation::EQ;
        }
      }
    }
    feasibility_->remove(interiorFace.cuts());
  }

  // Assign whether the input domain contains the origin
  BinaryVector origin(voronoi_.dimPoints(), 0);
  domainContainsOrigin_ = true;
  for (const auto& [a, b, rel] : domain_) {
    double value = dot(a, origin) + b;
    switch (rel) {
      case Relation::LT:
        if (value >= -params_.tolerance) {
          domainContainsOrigin_ = false;
          break;
        }
      case Relation::LE:
        if (value >= 0.0) {
          domainContainsOrigin_ = false;
          break;
        }
      case Relation::EQ:
        if (std::abs(value) > params_.tolerance) {
          domainContainsOrigin_ = false;
          break;
        }
      case Relation::GE:
        if (value <= -params_.tolerance) {
          domainContainsOrigin_ = false;
          break;
        }
      case Relation::GT:
        if (value <= 0.0) {
          domainContainsOrigin_ = false;
          break;
        }
      case Relation::RT:
        // Always true constraint, do nothing
        break;
      case Relation::RF:
        domainContainsOrigin_ = false;
        break;
      default:
        throw std::invalid_argument("Unknown relation in domain constraints.");
    }
    if (!domainContainsOrigin_) { break; }
  }
}

void Dynprog::initRootState() {
  // Initialize root state
  State root = State{Cone()};

  // Initialize root faces (all faces if no input domain, else filtered)
  root.faceIds.reserve(voronoi_.numFaces());
  if (domain_.empty()) {
    for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) {
      root.faceIds.push_back(faceId);
    }
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
  for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) {
    root.splits.push_back(Split{splitId});
  }

  // Initialize root lower bound
  evaluateLowerBound(root);

  // Register root state
  rootId_ = states_.size();
  regionToStateId_[root.region] = rootId_;
  states_.push_back(root);
}

State Dynprog::createState(State& parent, Split& split, Relation cutDir) {
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

void Dynprog::buildState(Index stateId) {
  stats_.numStatesBuilt++;

  // Temporary data <splitPos, splitId, dirId, childState> of new children. New
  // children are registered in bulk at the end of the method to avoid dangling
  // references of states/splits and reduce memory footprint in GREEDY mode.
  std::vector<std::tuple<Index, Index, Index, State>> newChildren;

  State& state = states_[stateId];
  std::vector<Split>& splits = state.splits;

  // ----- Build splits sequentially ----- //

  double bestSplitScore = std::numeric_limits<double>::infinity();
  Index bestSplitId = INVALID_INDEX;
  Index bestSplitLb = MAX_DEPTH;

  feasibility_->add(state.region.cuts());
  for (Index splitPos = 0; splitPos < splits.size(); ++splitPos) {
    Split& split = splits[splitPos];

    // Recover child indices from memoization or mark as not found
    for (Index i = 0; i < branchDirections_.size(); ++i) {
      Relation cutDir = branchDirections_[i];
      Cut childCut = Cut(split.id, cutDir);
      Cone childRegion = state.region.refine(childCut);
      auto itChild = regionToStateId_.find(childRegion);
      if (itChild != regionToStateId_.end()) {
        split.childIds[i] = itChild->second;
      } else {
        split.childIds[i] = INVALID_INDEX;
      }
    }

    // Check split validity
    evaluateValidity(state, split, true);

    // Continue if split is invalid
    if (!split.valid.value()) { continue; }

    // Recover child number of faces and lower bounds
    // Note: create missing child states and insert in temporary buffer
    Index bestChildLb = 0;
    std::vector<Index> childFaceCounts = {};
    for (Index dirId = 0; dirId < branchDirections_.size(); ++dirId) {
      Relation cutDir = branchDirections_[dirId];
      Index childId = split.childIds[dirId];
      if (childId != INVALID_INDEX) {
        bestChildLb = std::max(bestChildLb, safeAdd(states_[childId].lbHeight));
        childFaceCounts.push_back(states_[childId].faceIds.size());
      } else {
        State child = createState(state, split, cutDir);

        bestChildLb = std::max(bestChildLb, safeAdd(child.lbHeight));
        childFaceCounts.push_back(child.faceIds.size());

        newChildren.push_back(
            std::make_tuple(splitPos, split.id, dirId, std::move(child)));
      }
    }

    // Evaluate split score
    evaluateScore(state, split, childFaceCounts);

    // Update best split score, id, and lower bound
    if (split.score < bestSplitScore) {
      bestSplitScore = split.score;
      bestSplitId = split.id;
    }
    if (bestChildLb < bestSplitLb) { bestSplitLb = bestChildLb; }
  }
  feasibility_->remove(state.region.cuts());

  // ----- Insert new child states in memoization ----- //

  // Warning: from here, state and splits references are now dangling

  for (const auto& [splitPos, splitId, dirId, child] : newChildren) {
    // In GREEDY exploration mode, only create children for the best split
    if (params_.exploration == Exploration::GREEDY) {
      if (splitId != bestSplitId) { continue; }
    }

    // Register new child in memoization
    Index childId = states_.size();
    regionToStateId_[child.region] = childId;
    states_[stateId].splits[splitPos].childIds[dirId] = childId;
    states_.push_back(std::move(child));
  }

  // ----- Register building results in state ----- //

  // Register that the state has been built
  states_[stateId].isBuilt = true;

  // Remove invalid splits
  states_[stateId].splits.erase(
      std::remove_if(std::begin(states_[stateId].splits),
                     std::end(states_[stateId].splits),
                     [](const auto& split) { return !split.valid.value(); }),
      std::end(states_[stateId].splits));

  // Check if the state became a leaf (no valid splits)
  if (states_[stateId].splits.size() == 0) {
    states_[stateId].lbHeight = 0;
    states_[stateId].ubHeight = 0;
    states_[stateId].isClosed = true;
    stats_.numStatesLeafed++;
    return;
  }

  // Sort splits by ascending score
  std::stable_sort(std::begin(states_[stateId].splits),
                   std::end(states_[stateId].splits),
                   [](const Split& splitA, const Split& splitB) {
                     return splitA.score < splitB.score;
                   });

  // Update the state lower bound
  if (bestSplitLb != MAX_DEPTH) {
    states_[stateId].lbHeight =
        std::max(states_[stateId].lbHeight, bestSplitLb);
  }

  // Check if the state can be pruned with improved lower bound
  if (states_[stateId].depth() + states_[stateId].lbHeight >=
      states_[rootId_].ubHeight) {
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
  if (states_[stateId].depth() + states_[stateId].lbHeight >=
      states_[rootId_].ubHeight) {
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
  for (Index splitPos = 0; splitPos < states_[stateId].splits.size();
       ++splitPos) {
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

      // Update the split lower bound
      lbHeightSplit =
          std::max(lbHeightSplit, safeAdd(states_[childId].lbHeight));

      // Check if the split can be pruned
      if (states_[stateId].depth() + lbHeightSplit >= states_[0].ubHeight) {
        states_[stateId].splits[splitPos].isClosed = true;
        break;
      }

      // Process the child
      evaluateState(childId, kSplits);

      // Update the split upper bound
      ubHeightSplit =
          std::max(ubHeightSplit, safeAdd(states_[childId].ubHeight));

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
    if (lbHeightSplit == ubHeightSplit) {
      states_[stateId].splits[splitPos].isClosed = true;
    }
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
      static_cast<Index>(std::ceil(std::log(static_cast<double>(numFaces)) /
                                   std::log(static_cast<double>(base))));
}

void Dynprog::evaluateScore(State& state, Split& split,
                            const std::vector<Index>& childFaceCounts) {
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
      for (Index count : childFaceCounts) {
        maxCount = std::max(maxCount, count);
      }
      score = static_cast<double>(maxCount);
      break;
    }

    // Random scoring in [0,1)
    case SplitScoring::RANDOM: {
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

void Dynprog::evaluateValidity(State& state, Split& split,
                               bool externalStateCuts) {
  // Try to infer feasibility based on filtering rules
  if (params_.filterChecks) {
    split.valid = inferValidity(state, split);
    if (split.valid.has_value()) { return; }
  }

  // If filtering is not enabled or failed, perform full feasibility check
  Cut splitCut = Cut(split.id, Relation::EQ);
  if (!externalStateCuts) { feasibility_->add(state.region.cuts()); }
  feasibility_->add(splitCut);
  split.valid = feasibility_->check();
  feasibility_->remove(splitCut);
  if (!externalStateCuts) { feasibility_->remove(state.region.cuts()); }
}

bool Dynprog::checkChildFace(State& state, Split& split, Relation cutDir,
                             Index faceId, bool externalChildCuts) {
  // Try to infer whether the face is a child one based on filtering rules
  if (params_.filterChecks) {
    std::optional<bool> isInferredChildFace =
        inferChildFace(state, split, cutDir, faceId);
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

std::optional<bool> Dynprog::inferValidity(State& state, Split& split) {
  // A split is certified valid if:
  // - at least two children are non-empty (separation occurs)
  // A split is certified invalid if:
  // - at least one child has identical faces (no separation)

  Index stateNumFaces = state.faceIds.size();
  Index childNonempty = 0;
  Index childIdentical = 0;

  for (Index i = 0; i < branchDirections_.size(); ++i) {
    Relation cutDir = branchDirections_[i];
    Index childId = split.childIds[i];
    Index childNumFaces = 0;

    // Retrieve child number of faces either:
    // - exactly from memoization
    // - as an under-approximation by inferring face intersections
    if (childId != INVALID_INDEX) {
      childNumFaces = states_[childId].faceIds.size();
    } else {
      childNumFaces = 0;
      for (Index faceId : state.faceIds) {
        std::optional<bool> isEmpty =
            inferChildFace(state, split, cutDir, faceId);
        if (isEmpty.has_value() && isEmpty.value()) { childNumFaces += 1; }
      }
    }

    if (childNumFaces > 0) { childNonempty += 1; }
    if (childNumFaces == stateNumFaces) { childIdentical += 1; }

    if (childNonempty >= 2) { return true; }
    if (childIdentical >= 1) { return false; }
  }

  // Fallback: no inference possible
  return std::nullopt;
}

std::optional<bool> Dynprog::inferChildFace(State& state, Split& split,
                                            Relation cutDir, Index faceId) {
  Relation position = positions_[faceId][split.id];
  const Cone& stateCone = state.region;

  if (position == Relation::LT) {
    if (cutDir == Relation::GT) {
      return false;
    } else if (cutDir == Relation::LT) {
      if (stateCone.isOpen()) { return true; }
    } else if (cutDir == Relation::EQ) {
      if (stateCone.containsOrigin() && domainContainsOrigin_) { return true; }
    }
  } else if (position == Relation::GT) {
    if (cutDir == Relation::LT) {
      return false;
    } else if (cutDir == Relation::GT) {
      if (stateCone.isOpen()) { return true; }
    } else if (cutDir == Relation::EQ) {
      if (stateCone.containsOrigin() && domainContainsOrigin_) { return true; }
    }
  } else if (position == Relation::EQ) {
    ///
  }

  if (childContainsFaceCenter(state, split, cutDir, faceId)) { return true; }

  return std::nullopt;
}

bool Dynprog::childContainsFaceCenter(State& state, Split& split,
                                      Relation cutDir, Index faceId) {
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
      for (Index k = 1; k <= voronoi_.numSplits(); ++k) {
        iterationRange.push_back(k);
      }
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
  double elapsed =
      std::chrono::duration<double>(Clock::now() - startTime_).count();
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
  if (elapsedTime(checkTime_) < params_.logInterval && message.empty()) {
    return;
  }

  checkTime_ = Clock::now();

  std::ostream& out = *(params_.outputStream);
  out << "  ";
  out << std::setw(12) << std::fixed << std::setprecision(2)
      << elapsedTime(startTime_);
  out << std::setw(12) << stats_.numIters;
  out << std::setw(12) << states_[rootId_].lbHeight;
  out << std::setw(12)
      << (states_[rootId_].ubHeight == MAX_DEPTH
              ? "inf"
              : std::to_string(states_[rootId_].ubHeight));
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
  out << "  time: " << std::fixed << std::setprecision(4) << stats_.runTime
      << "\n";
}

}  // namespace treeco
