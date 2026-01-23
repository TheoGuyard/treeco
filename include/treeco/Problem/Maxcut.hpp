/**
 * @file Maxcut.hpp
 * @brief Maximum-weighted cut problem with fixed number of vertices.
 */

#ifndef TREECO_PROBLEM_MAXCUT_HPP
#define TREECO_PROBLEM_MAXCUT_HPP

#include <algorithm>
#include <utility>

#include "treeco/Problem.hpp"

namespace treeco {

/**
 * @brief Maximum Cut Problem on a complete undirected graph.
 * 
 * Given d vertices, the problem is to partition vertices into two sets S and T
 * such that the total weight of edges crossing the partition is maximized.
 * 
 * Formulation as max <c, x>:
 *   - Variables: x_ij for edges (i,j) with i < j, so d*(d-1)/2 edge variables
 *   - x_ij = 1 if edge (i,j) is cut (endpoints in different partitions)
 *   - c_ij = weight of edge (i,j)
 *   - Feasible set X = all binary vectors corresponding to valid cuts
 * 
 * Edge indexing: edge (i,j) with i < j has index i*d - i*(i+1)/2 + (j-i-1)
 */
class Maxcut : public Problem {
public:
    /**
     * @brief Construct a MaxCut problem.
     * @param numVertices Number of vertices in the complete graph
     */
    explicit Maxcut(Index numVertices);

    /**
     * @brief Get the number of vertices.
     * @return Number of vertices in the graph.
     */
    Index numVertices() const { return numVertices_; }

    /**
     * @brief Get all feasible solutions (all valid cuts).
     * @return Vector of all feasible solutions, represented as x = {x_ij} for
     * all edges (i,j) with i < j.
     */
    std::vector<BinaryVector> getFeasibleSet() const override;

    /**
     * @brief Get the cost vector domain (positive entires).
     * @return Domain of the cost vector C = {c | c_ij > 0 for all edges (i,j)}
     */
    Domain getCostDomain() const;

    using Problem::sampleCost;

    /**
     * @brief Sample a random cost vector using a provided seed.
     * @param seed Seed for the random number generator
     * @return Cost vector with entries uniformly distributed in (0, 1]
     */
    RealVector sampleCost(Index seed) const;

    /**
     * @brief Convert edge indices to variable index.
     * @param i First vertex (i < j)
     * @param j Second vertex
     * @return Variable index for edge (i,j)
     */
    Index edgeToIndex(Index i, Index j) const;

    /**
     * @brief Convert variable index to edge indices.
     * @param idx Variable index
     * @return Pair (i, j) with i < j
     */
    std::pair<Index, Index> indexToEdge(Index idx) const;

private:
    Index numVertices_; // Number of vertices
};

} // namespace treeco

#endif // TREECO_PROBLEM_MAXCUT_HPP
