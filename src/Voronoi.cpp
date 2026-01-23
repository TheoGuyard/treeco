#include "treeco/Voronoi.hpp"

namespace treeco {

Voronoi::Voronoi(const std::vector<SimplexVector>& points) : points_(points) {
    if (points_.empty()) {
        throw std::invalid_argument("Point list cannot be empty");
    }
    if (points_.size() < 2) {
        throw std::invalid_argument("Point list contains less than 2 points");
    }
    for (const auto& p : points_) {
        if (p.size() != points_[0].size()) {
            throw std::invalid_argument("Points must have the same dimension");
        }
    }
}

void Voronoi::build(const VoronoiParams& params) {

    clear();

    params_ = params;
    startTime_ = Clock::now();
    checkTime_ = startTime_;

    Adjacency adj(points_);
    
    logHeader();

    // Prepare memory allocation
    faces_.resize(numPoints());
    edges_.reserve(numPoints() * (numPoints() - 1));
    splits_.reserve(numPoints() * (numPoints() - 1));

    // Main construction loop
    std::map<TernaryVector, Index> splitToIndex  = {};
    for (Index i = 0; i < numPoints(); ++i) {

        faces_[i].pointId = i;

        for (Index j = i + 1; j < numPoints(); ++j) {
        
            logProgress(adj, i);

#ifdef TREECO_BUILD_PYTHON
            // Propagate keyboard interrupts from python
            py::gil_scoped_acquire acquire;
            if (PyErr_CheckSignals() != 0) { throw py::error_already_set(); }
#endif

            if (adj.check(i, j)) {
                
                // Compute split
                TernaryVector split = bisector(points_[i], points_[j]);

                // Check split side
                int splitSide = dot(split, points_[i]);

                // Check if split already exists or insert it
                auto [it, inserted] = splitToIndex.try_emplace(split, splitToIndex.size());
                if (inserted) { splits_.push_back(std::move(split)); }
                int splitId = it->second;

                // Check on which side each face lies
                if (splitSide > 0) {
                    faces_[i].cone.addCut(splitId, Relation::GE);
                    faces_[j].cone.addCut(splitId, Relation::LE);
                    edges_.emplace_back(splitId, j, i);
                } else {
                    faces_[i].cone.addCut(splitId, Relation::LE);
                    faces_[j].cone.addCut(splitId, Relation::GE);
                    edges_.emplace_back(splitId, i, j);
                }
            }
        }
    }

    logProgress(adj, numPoints(), "done");

    // Shrink back to fit actual size
    faces_.shrink_to_fit();
    edges_.shrink_to_fit();
    splits_.shrink_to_fit();

    // Finalize stats
    stats_.buildTime = std::chrono::duration<double>(Clock::now() - startTime_).count();
    stats_.lpSolved = adj.getNumSolve();
    stats_.isBuilt = true;

    logFooter();
}

void Voronoi::clear() {
    splits_.clear();
    faces_.clear();
    edges_.clear();
    stats_ = VoronoiStats();
}

void Voronoi::logHeader() const {
    if (!params_.verbose) { return; }
    std::ostream& out = *(params_.outputStream);
    out << "Constructing Voronoi diagram...\n";
    out << "  dim points: " << dimPoints() << "\n";
    out << "  num points: " << numPoints() << "\n";
    out << std::string("  ") + std::string(5 * 12, '-') << "\n";
    out << "  ";
    out << std::setw(12) << "time";
    out << std::setw(12) << "cells";
    out << std::setw(12) << "edges";
    out << std::setw(12) << "splits";
    out << std::setw(12) << "lp opt";
    out <<  "\n";
}

void Voronoi::logProgress(const Adjacency& adj, Index i, const std::string& message) {
    if (!params_.verbose) { return; }
    if (elapsedTime(checkTime_) < params_.logInterval && message.empty()) { return; }
    
    checkTime_ = Clock::now();

    std::ostream& out = *(params_.outputStream);
    out << "  ";
    out << std::setw(12) << std::fixed << std::setprecision(2) << elapsedTime(startTime_);
    out << std::setw(12) << i;
    out << std::setw(12) << edges_.size();
    out << std::setw(12) << splits_.size();
    out << std::setw(12) << adj.getNumSolve();
    out << "  " << message;
    out << "\n";                        
}

void Voronoi::logFooter() const {
    if (!params_.verbose) { return; }
    std::ostream& out = *(params_.outputStream);     
    out << std::string("  ") + std::string(5 * 12, '-') << "\n";   
    out << "  time: " << std::fixed << std::setprecision(4) << stats_.buildTime << "\n";    
}

std::ostream& operator<<(std::ostream& oss, const Voronoi& voronoi) {
    oss << "Voronoi diagram\n";
    oss << "  dimension : " << voronoi.dimPoints() << "\n";
    oss << "  num faces : " << voronoi.numFaces() << "\n";
    oss << "  num edges : " << voronoi.numEdges() << "\n";
    oss << "  num splits: " << voronoi.numSplits() << "\n";
    return oss;
}

} // namespace treeco
