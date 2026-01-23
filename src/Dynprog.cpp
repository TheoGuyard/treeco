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
    oss << "  lp solved      : " << stats.numFeasibilityChecks << "\n";
    if (stats.optimalDepth == MAX_DEPTH) {
        oss << "  optimal depth  : inf (not found)";
    } else {
        oss << "  optimal depth  : " << stats.optimalDepth << "";
    }
    return oss;
}


Dynprog::Dynprog(
    const Voronoi& voronoi,
    const Domain& domain
) : voronoi_(voronoi),
    domain_(domain)
{
    if (!voronoi_.isBuilt()) {
        throw std::invalid_argument("Voronoi diagram must be built before using Dynprog");
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

    // Initialize random number generator for split sampling
    rng_.seed(params_.randomSeed);
        
    // Initialize feasibility checker with splits and input domain specs
    feasibility_ = std::make_unique<Feasibility>(
        voronoi_.splits(),
        domain_,
        params_.tolerance
    );
    
    states_.clear();
    regionToStateId_.clear();
    initPositions();    
    initRootState();
    
    // ------------------------------------------
    // Iterative exploration loop
    // ------------------------------------------
    
    logProgress("initialization");

    for (Index kSplits : getIterationRange()) {
        stats_.numIters = kSplits;
        
        evaluateState(rootId_, kSplits);
        updateStatus();

        if (status_ == DynprogStatus::OPTIMAL) {
            logProgress("optimal");
            break;
        } else if (checkTimeLimit()) {
            logProgress("time limit");
            break;
        }
        
        assert(status_ == DynprogStatus::SUBOPTIMAL);
        logProgress("subopt (" + std::to_string(kSplits) + " splits)");
    }
        
    // ------------------------------------------
    // Finalization
    // ------------------------------------------

    // Fill remaining statistics fields
    stats_.runTime = std::chrono::duration<double>(Clock::now() - startTime_).count();
    stats_.numStates = states_.size();
    stats_.numFeasibilityChecks = feasibility_->numSolve();
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

    // Precompute all positions if requested
    if (params_.positioning == Positioning::PRECOMPUTE) {
        for (Index faceId = 0; faceId < voronoi_.numFaces(); ++faceId) {

#ifdef TREECO_BUILD_PYTHON
            // Propagate keyboard interrupts from python
            py::gil_scoped_acquire acquire;
            if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
#endif

            const Cone& interiorFace = voronoi_.face(faceId).cone.interior();
            feasibility_->add(interiorFace.cuts());
            for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) {
                positions_[faceId][splitId] = getPosition(faceId, splitId, true);
            }
            feasibility_->remove(interiorFace.cuts());
        }
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
            if (feasibility_->check()) {
                root.faceIds.push_back(faceId);
            }
            feasibility_->remove(face.cone.cuts());
        }
    }
    root.faceIds.shrink_to_fit();

    // Initialize root splits (all unbuild splits)
    root.splitIds.reserve(voronoi_.numSplits());
    for (Index splitId = 0; splitId < voronoi_.numSplits(); ++splitId) {
        root.splitIds.push_back(splitId);
        root.splits[splitId] = Split();
    }
    root.splitIds.shrink_to_fit();
    
    // Initialize root lower bound
    evaluateLowerBound(root);

    // Register root state
    rootId_ = states_.size();
    regionToStateId_[root.region] = rootId_;
    states_.push_back(root);
}

State Dynprog::createState(Index parentId, Index splitId, Relation cutDir) {

    const State& parent = states_[parentId];

    Cut cut = Cut(splitId, cutDir);
    Cone region = parent.region.refine(cut);
    State child = State{region};

    // Initialize child faces
    child.faceIds.reserve(parent.faceIds.size());
    feasibility_->add(child.region.cuts());
    for (Index faceId : parent.faceIds) {

#ifdef TREECO_BUILD_PYTHON
        // Propagate keyboard interrupts from python
        py::gil_scoped_acquire acquire;
        if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
#endif

        bool intersect = checkChildFace(parentId, splitId, cutDir, faceId, true);

        if (intersect) {
            child.faceIds.push_back(faceId);
        }
    }
    feasibility_->remove(child.region.cuts());
    child.faceIds.shrink_to_fit();

    // Initialize child splits
    child.splitIds.reserve(std::max(parent.splitIds.size() - 1, Index(0)));
    for (Index parentSplitId : parent.splitIds) {                
        if (parentSplitId == splitId) { continue; }
        child.splitIds.push_back(parentSplitId);
        child.splits[parentSplitId] = Split();
    }
    child.splitIds.shrink_to_fit();

    // Initialize child lower bound
    evaluateLowerBound(child);

    // Test if new child is a leaf
    if (child.isLeaf()) {
        child.lbHeight = 0;
        child.ubHeight = 0;
        child.isClosed = true;
    }
    
    return child;
}

void Dynprog::buildState(Index stateId, Index kSplits) {

    stats_.numStatesBuilt++;

    Index stateFaceCount = states_[stateId].faceIds.size();

    // ----- Initialize building process ----- //

    // Temporary data of new children: newChildren[splitId][cutDir] = State
    // Note: All new children are registered in bulk at the end of the method
    std::map<Index, std::map<Relation, State>> newChildren;

    // Retrieve candidate splits not checked for validity yet
    std::vector<Index> validSplitsIds;
    std::vector<Index> candidateSplitsIds;
    validSplitsIds.reserve(states_[stateId].splitIds.size());
    candidateSplitsIds.reserve(states_[stateId].splitIds.size() - states_[stateId].numSplitsBuilt);
    for (Index splitId : states_[stateId].splitIds) {
        if (!states_[stateId].splits[splitId].valid.has_value()) {
            candidateSplitsIds.push_back(splitId);
        } else {
            assert(states_[stateId].splits[splitId].valid.value());
            validSplitsIds.push_back(splitId);
        }
    }

    // Shuffle candidate splits if sampling is activated
    if (params_.splitSelection == SplitSelection::SAMPLING) {
        std::shuffle(
            candidateSplitsIds.begin(),
            candidateSplitsIds.end(),
            rng_
        );
    }


    // ----- Build splits sequentially ----- //

    Index numCandidateChecked = 0;
    Index bestSplitLb = MAX_DEPTH;
    feasibility_->add(states_[stateId].region.cuts());
    for (Index splitId : candidateSplitsIds) {

        // For SplitSelection::SAMPLING, if kSplits have already been built,
        // simply add remaining candidates to the valid splits pool
        if (params_.splitSelection == SplitSelection::SAMPLING) {
            if (states_[stateId].numSplitsBuilt >= kSplits) {
                validSplitsIds.push_back(splitId);
                continue;
            }
        }

        numCandidateChecked++;

        Split& split = states_[stateId].splits[splitId];

        // Recover child indices from memoization or mark as not found
        for (Relation cutDir : branchDirections_) {
            Cut childCut = Cut(splitId, cutDir);
            Cone childRegion = states_[stateId].region.refine(childCut);
            auto itChild = regionToStateId_.find(childRegion);
            if (itChild != regionToStateId_.end()) {
                split.childIds[cutDir] = itChild->second;
            } else {
                split.childIds[cutDir] = INVALID_INDEX;
            }
        }

        // Check split validity
        split.valid = checkSplitValidity(stateId, splitId, true);

        assert(split.valid.has_value());

        // Remove split data if not valid and continue
        if (!split.valid.value()) {
            states_[stateId].splits.erase(splitId);
            continue;
        }

        // Register valid split
        validSplitsIds.push_back(splitId);

        // Recover child number of faces or create new child states and update
        // the best split lower bound
        for (const auto& [cutDir, childId] : split.childIds) {
            if (childId != INVALID_INDEX) {
                split.childNumFaces[cutDir] = states_[childId].faceIds.size();
                bestSplitLb = std::min(
                    bestSplitLb,
                    safeAdd(states_[childId].lbHeight)
                );
            } else {
                State child = createState(stateId, splitId, cutDir);
                split.childNumFaces[cutDir] = child.faceIds.size();
                bestSplitLb = std::min(
                    bestSplitLb,
                    safeAdd(child.lbHeight)
                );

                // Append to temporary buffer
                newChildren[splitId][cutDir] = std::move(child);
            }
        }

        // Evaluate split score
        evaluateSplitScore(stateId, splitId);

        // Register split built
        states_[stateId].numSplitsBuilt++;
    }
    feasibility_->remove(states_[stateId].region.cuts());


    // ----- Register building results in state ----- //

    // Set back valid splits ids
    states_[stateId].splitIds = std::move(validSplitsIds);

    // Sort splits id by score
    const std::map<Index, Split>& stateSplits = states_[stateId].splits;
    std::sort(
        states_[stateId].splitIds.begin(),
        states_[stateId].splitIds.end(),
        [&stateSplits](Index a, Index b) {
            return stateSplits.at(a).score < stateSplits.at(b).score;
        }
    );

    // Mark state as built if all candidates have been processed
    if (numCandidateChecked >= candidateSplitsIds.size()) {
        states_[stateId].isBuilt = true;

        // Update the state lower bound
        if (params_.lowerBounding == LowerBounding::BACKTRACK) {
            if (bestSplitLb != MAX_DEPTH) {
                states_[stateId].lbHeight = std::max(
                    states_[stateId].lbHeight,
                    bestSplitLb
                );
            }
        }

        // Check if the state became a leaf (no valid splits)
        if (states_[stateId].numSplitsBuilt == 0) {
            states_[stateId].lbHeight = 0;
            states_[stateId].ubHeight = 0;
            states_[stateId].isClosed = true;
            stats_.numStatesLeafed++;
        }

        // Check if the state can be pruned (with improved lower bound)
        else if (states_[stateId].depth() + states_[stateId].lbHeight >= states_[rootId_].ubHeight) {
            states_[stateId].lbHeight = MAX_DEPTH;
            states_[stateId].ubHeight = MAX_DEPTH;
            states_[stateId].isClosed = true;
            stats_.numStatesPruned++;
        }
    }


    // ----- Insert new child states in memoization ----- //

    for (const auto& [splitId, splitChildren] : newChildren) {

        // In GREEDY exploration mode, only create children for the first split
        if (params_.exploration == Exploration::GREEDY) {
            if (splitId != states_[stateId].splitIds[0]) {
                continue;
            }
        }

        // Register new child in memoization
        for (const auto& [cutDir, child] : splitChildren) {
            Index childId = states_.size();
            regionToStateId_[child.region] = childId;
            states_[stateId].splits[splitId].childIds[cutDir] = childId;
            states_.push_back(std::move(child));
        }
    }
}

void Dynprog::evaluateState(Index stateId, Index kSplits) {

#ifdef TREECO_BUILD_PYTHON
    // Propagate keyboard interrupts from python
    py::gil_scoped_acquire acquire;
    if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
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
    
    // Check if the state splits need to be built
    if (!states_[stateId].isBuilt && states_[stateId].numSplitsBuilt < kSplits) {
        buildState(stateId, kSplits);        
        // Check if state has been closed during building
        if (states_[stateId].isClosed) { return; }
    }

    // -------------------- State processing -------------------- //
    
    Index numSplitsTested = 0;
    Index bestSplitLb = MAX_DEPTH;
    for (Index splitId : states_[stateId].splitIds) {

        // Break if enough splits have been tested
        if (numSplitsTested >= kSplits) { break; }

        // Skip if the split is closed (i.e., pruned, or optimal)
        if (states_[stateId].splits[splitId].isClosed) { continue; }

        numSplitsTested++;

        // Process each children
        Index lbHeightSplit = 0;
        Index ubHeightSplit = 0;
        for (const auto& [cutDir, childId] : states_[stateId].splits[splitId].childIds) {

            // Update the split lower bound
            lbHeightSplit = std::max(lbHeightSplit, safeAdd(states_[childId].lbHeight));

            // Check if the split can be pruned
            if (states_[stateId].depth() + lbHeightSplit >= states_[0].ubHeight) {
                states_[stateId].splits[splitId].isClosed = true;
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
        if (states_[stateId].splits[splitId].isClosed) { continue; }

        // Update best split lower bound
        bestSplitLb = std::min(bestSplitLb, lbHeightSplit);

        // Update the best state upper bound and split index                
        if (ubHeightSplit < states_[stateId].ubHeight) {
            states_[stateId].splitId = splitId;
            states_[stateId].ubHeight = ubHeightSplit;
        }
        
        // Check if the split is optimal by matching bounds
        if (lbHeightSplit == ubHeightSplit) {
            states_[stateId].splits[splitId].isClosed = true;
        }
    }

    // Update state lower bound if backtracking
    if (params_.lowerBounding == LowerBounding::BACKTRACK) {
        if (bestSplitLb != MAX_DEPTH) {
            states_[stateId].lbHeight = std::max(
                states_[stateId].lbHeight,
                bestSplitLb
            );
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
    if (kSplits >= states_[stateId].splitIds.size()) {
        states_[stateId].lbHeight = states_[stateId].ubHeight;
        states_[stateId].isClosed = true;
        stats_.numStatesClosed++;
        return;
    }
}

void Dynprog::evaluateLowerBound(State& state) {    
    Index numFaces = state.faceIds.size();

    if (numFaces <= 1) {
        state.lbHeight = 0;
        return;
    }
    
    Index base = branchDirections_.size();

    state.lbHeight = static_cast<Index>(
        std::ceil(
            std::log(static_cast<double>(numFaces)) / 
            std::log(static_cast<double>(base))
        )
    );
}

void Dynprog::evaluateSplitScore(Index stateId, Index splitId) {

    State& state = states_[stateId];
    Split& split = state.splits[splitId];
    
    Index stateFaceCount = state.faceIds.size();
    Index numChildren = split.childIds.size();
    double score;
    
    switch (params_.splitScoring) {

        // Variance from equal split: lower is better
        case SplitScoring::VARIANCE: {
            double mean = static_cast<double>(stateFaceCount) / numChildren;
            double svar = 0.0;
            for (const auto& [_, count] : split.childNumFaces) {
                double diff = static_cast<double>(count) - mean;
                svar += diff * diff;
            }
            score = svar / numChildren;
            break;
        }
        
        // Information gain: higher reduction is better (set negative)
        case SplitScoring::ENTROPY: {
            Index totalChildFaces = 0;
            for (const auto& [_, count] : split.childNumFaces) {
                totalChildFaces += count;
            }
            if (totalChildFaces == 0) {
                score = std::numeric_limits<double>::infinity();
                break;
            }
            double stateEntropy = std::log2(static_cast<double>(stateFaceCount));
            double childEntropy = 0.0;
            for (const auto& [_, count] : split.childNumFaces) {
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
            for (const auto& [_, count] : split.childNumFaces) {
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

void Dynprog::updateStatus() {
    if (states_[rootId_].isClosed) {
        status_ = DynprogStatus::OPTIMAL;
    } else if (states_[rootId_].ubHeight < MAX_DEPTH) {
        status_ = DynprogStatus::SUBOPTIMAL;
    } else {
        status_ = DynprogStatus::INVALID;
    }
}

Relation Dynprog::getPosition(Index faceId, Index splitId, bool externalInteriorFaceCuts) {
    
    // Check cache first
    if (positions_[faceId][splitId] != Relation::RF) {
        return positions_[faceId][splitId];
    }

    const Cone& interiorCone = voronoi_.face(faceId).cone.interior();
    
    // Check if the split if a face boundary
    for (const auto& [sid, dir] : interiorCone.cuts()) {
        if (sid == splitId) {
            positions_[faceId][splitId] = dir;
            return dir;
        }
    }

    // Determine position with feasibility check if not boundary
    // - int(face) ∩ s< = ∅ → face ⊆ s>=
    // - int(face) ∩ s> = ∅ → face ⊆ s<=
    // - otherwise → int(face) intersects s< and s>

    if (!externalInteriorFaceCuts) { feasibility_->add(interiorCone.cuts()); }

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

    if (!externalInteriorFaceCuts) { feasibility_->remove(interiorCone.cuts()); }

    return positions_[faceId][splitId];
}

bool Dynprog::checkSplitValidity(Index stateId, Index splitId, bool externalStateCuts) {
    
    // Try to infer feasibility based on filtering rules
    if (params_.filterChecks) {
        std::optional<bool> isInferredValid = inferSplitValidity(stateId, splitId);
        if (isInferredValid.has_value()) {
            return isInferredValid.value();
        }
    }

    // If filtering is not enabled or failed, perform full feasibility check
    Cut splitCut = Cut(splitId, Relation::EQ);
    if (!externalStateCuts) { feasibility_->add(states_[stateId].region.cuts()); }
    feasibility_->add(splitCut);
    bool isValid = feasibility_->check();
    feasibility_->remove(splitCut);
    if (!externalStateCuts) { feasibility_->remove(states_[stateId].region.cuts()); }
    
    return isValid;
}

bool Dynprog::checkChildFace(Index stateId, Index splitId, Relation cutDir, Index faceId, bool externalChildCuts) {

    // Try to infer whether the face is a child one based on filtering rules
    if (params_.filterChecks) {
        std::optional<bool> isInferredChildFace = inferChildFace(stateId, splitId, cutDir, faceId);
        if (isInferredChildFace.has_value()) { return isInferredChildFace.value(); }
    }
    
    // If filtering is not enabled or failed, perform full feasibility check
    if (!externalChildCuts) {
        feasibility_->add(states_[stateId].region.cuts());
        feasibility_->add(Cut(splitId, cutDir));
    }
    feasibility_->add(voronoi_.face(faceId).cone.cuts());
    bool isChildFace = feasibility_->check();
    feasibility_->remove(voronoi_.face(faceId).cone.cuts());
    if (!externalChildCuts) {
        feasibility_->remove(states_[stateId].region.cuts());
        feasibility_->remove(Cut(splitId, cutDir));
    }

    return isChildFace;
}

std::optional<bool> Dynprog::inferSplitValidity(Index stateId, Index splitId) {

    // A split is certified valid if:
    // - at least two children are non-empty (separation occurs)
    // A split is certified invalid if:
    // - at least one child has identical faces (no separation) 

    Index stateNumFaces = states_[stateId].faceIds.size();
    Index childNonempty = 0;
    Index childIdentical = 0;

    for (const auto& [cutDir, childId] : states_[stateId].splits.at(splitId).childIds) {

        Index childNumFaces = 0;
        
        // Retrieve child number of faces either:
        // - exactly from memoization
        // - as an under-approximation by inferring face intersections
        if (childId != INVALID_INDEX) {
            childNumFaces = states_[childId].faceIds.size();
        } else {
            childNumFaces = 0;
            for (Index faceId : states_[stateId].faceIds) {
                std::optional<bool> isEmpty = inferChildFace(stateId, splitId, cutDir, faceId);
                if (isEmpty.has_value() && isEmpty.value()) {
                    childNumFaces += 1;
                }
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

std::optional<bool> Dynprog::inferChildFace(Index stateId, Index splitId, Relation cutDir, Index faceId) {

    Relation position = getPosition(faceId, splitId);
    const Cone& stateCone = states_[stateId].region;

    if (stateCone.isOpen()) {
        if (position == Relation::LT) {
            return (cutDir == Relation::LT || cutDir == Relation::EQ);
        } else if (position == Relation::GT) {
            return (cutDir == Relation::GT || cutDir == Relation::EQ);
        } else if (position == Relation::EQ) {
            return true;
        }
    }

    if (position == Relation::LT) {
        if (cutDir == Relation::LT) {
            if (stateCone.isOpen()) {
                return true;
            }
        } else if (cutDir == Relation::EQ) {
            if (stateCone.isOpen()) {
                return true;
            } else if (stateCone.containsOrigin() && domainContainsOrigin_) {
                return true;
            }
        } else if (cutDir == Relation::GT) {
            return false;
        } else {
            throw std::runtime_error("Invalid cut relation.");
        }
    } else if (position == Relation::GT) {
        if (cutDir == Relation::LT) {
            return false;
        } else if (cutDir == Relation::EQ) {
            if (stateCone.isOpen()) {
                return true;
            } else if (stateCone.containsOrigin() && domainContainsOrigin_) {
                return true;
            }
        } else if (cutDir == Relation::GT) {
            if (stateCone.isOpen()) {
                return true;
            }
        }  else {
            throw std::runtime_error("Invalid cut relation.");
        }
    } else if (position == Relation::EQ) {
        if (childContainsFaceCenter(stateId, splitId, cutDir, faceId)) {
            return true;
        }
    } else {
        throw std::runtime_error("Invalid position.");
    }

    // Fallback: no inference possible
    return std::nullopt;
}

bool Dynprog::childContainsFaceCenter(Index stateId, Index splitId, Relation cutDir, Index faceId) {

    const Cone& childCone = states_[stateId].region.refine(Cut(splitId, cutDir));
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

void Dynprog::logFooter() const {
    if (!params_.verbose) return;
    std::ostream& out = *(params_.outputStream);
    out << std::string("  ") + std::string(9 * 12, '-') << "\n";
    out << "  time: " << std::fixed << std::setprecision(4) << stats_.runTime << "\n";
}

} // namespace treeco
