import numpy as np
from treeco import Tsp, LDTree, SplitScoring

# Create a TSP problem with 4 cities
problem = Tsp(num_cities=4)

# Extract all possible feasible points
points = problem.get_feasible_set()

# Build a LDTree policy for this problem
policy = LDTree(points)
policy.build(verbose=True)

policy = LDTree(points)
policy.build(verbose=True, split_scoring=SplitScoring.minmax)


# Print the tree structure
print("\n\n")
policy.pprint()
print("\n\n")

# Perform queries for various cost vectors
for i in range(10):
    cost = problem.sample_cost()
    sols = policy.query(cost)
    vals = np.array([np.dot(cost, sol) for sol in sols])
    print(f"query {i}")
    print(f"  cost : {cost}")
    print(f"  sols : {sols}")
    print(f"  vals : {vals}")
