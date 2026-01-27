import numpy as np
from treeco import LDTree
from treeco.problem import Tsp

# Create a TSP problem with 4 cities
problem = Tsp(num_cities=4)

# Extract all possible feasible points
points = problem.get_feasible_set()

# Extract the cost vector domain
domain = problem.get_cost_domain()

# Build a LDTree policy for this problem
policy = LDTree(points, domain)
policy.build(verbose=True)

policy = LDTree(points)
policy.build(verbose=True)


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
