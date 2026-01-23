#include <cstring>
#include <limits>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <treeco/LDTree.hpp>
#include <treeco/Problem/Tsp.hpp>
#include <treeco/Types.hpp>

using namespace treeco;

int main(int argc, char* argv[]) {
    
    // Create a TSP problem with 4 cities
    Tsp problem(4);

    // Extract all possible feasible points
    std::vector<BinaryVector> points = problem.getFeasibleSet();

    // Build a LDTree policy for this problem
    LDTree policy(points);
    policy.build(true);

    // Print the tree structure
    std::cout << "\n\n";
    policy.pprint();
    std::cout << "\n\n";

    // Perform queries for various cost vectors
    for (size_t i = 0; i < 10; ++i) {
        RealVector cost = problem.sampleCost();
        std::vector<BinaryVector> sols = policy.query(cost);

        std::cout << "query " << (i + 1) << "\n";
        std::cout << "  cost : ";
        printVector(cost, &std::cout);
        std::cout << "\n";
        std::cout << "  sols : {";
        for (const BinaryVector& sol : sols) {
            printVector(sol, &std::cout);
            if (&sol != &sols.back()) std::cout << ",";
        }
        std::cout << "}\n";
        std::cout << "  vals : {";
        for( const BinaryVector& sol : sols) {
            double val = dot(cost, sol);
            std::cout << val;
            if (&sol != &sols.back()) std::cout << ",";
        }
        std::cout << "}\n";
    }
    
    return 0;
}
