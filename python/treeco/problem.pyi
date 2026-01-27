from typing import List, Tuple, overload

from treeco.geometry import Domain
from treeco.types import BinaryVector, RealVector

class Problem:
    """
    Abstract base class for optimization problems.
    """

    @property
    def dimension(self) -> int:
        """Dimension of the problem."""
        ...

    def get_feasible_set(self) -> List[BinaryVector]:
        """Get all feasible solutions as a list of binary vectors."""
        ...

    def get_cost_domain(self) -> Domain:
        """Get the domain constraints for the cost vector."""
        ...

    @overload
    def sample_cost(self) -> RealVector:
        """Sample a random cost vector with entries uniformly in [0, 1]."""
        ...

    @overload
    def sample_cost(self, seed: int) -> RealVector:
        """Sample a random cost vector using a provided seed."""
        ...

class Explicit(Problem):
    """
    Optimization problem with an explicitly provided feasible set.
    """

    def __init__(self, feasible_set: List[BinaryVector]) -> None:
        """
        Construct from an explicit feasible set.

        Args:
            feasible_set: List of all feasible binary solutions
        """
        ...

class Knapsack(Problem):
    """
    0/1 Knapsack Problem.
    """

    @overload
    def __init__(self, num_items: int) -> None:
        """
        Construct a Knapsack problem with given number n of items, weights
        {1,...,n}, and capacity C = n.

        Args:
            num_items: Number of items
        """
        ...

    @overload
    def __init__(self, weights: List[float], capacity: float) -> None:
        """
        Construct a Knapsack problem with given item weights and capacity.

        Args:
            weights: List of item weights
            capacity: Maximum weight capacity of the knapsack
        """
        ...

    @property
    def num_items(self) -> int:
        """Get the number of items."""
        ...

    @property
    def weights(self) -> RealVector:
        """Get the item weights."""
        ...

    @property
    def capacity(self) -> float:
        """Get the knapsack capacity."""
        ...

class Maxcut(Problem):
    """
    Maximum Cut Problem on a complete undirected graph.
    """

    def __init__(self, num_vertices: int) -> None:
        """
        Construct a MaxCut problem.

        Args:
            num_vertices: Number of vertices in the complete graph
        """
        ...

    @property
    def num_vertices(self) -> int:
        """Get the number of vertices."""
        ...

    def edge_to_index(self, i: int, j: int) -> int:
        """
        Convert edge indices to variable index.

        Args:
            i: First vertex (i < j)
            j: Second vertex

        Returns:
            Variable index for edge (i,j)
        """
        ...

    def index_to_edge(self, idx: int) -> Tuple[int, int]:
        """
        Convert variable index to edge indices.

        Args:
            idx: Variable index

        Returns:
            Tuple (i, j) with i < j
        """
        ...

class Tsp(Problem):
    """
    Traveling Salesman Problem on a complete undirected graph.
    """

    def __init__(self, num_cities: int) -> None:
        """
        Construct a TSP instance.

        Args:
            num_cities: Number of cities in the complete graph
        """
        ...

    @property
    def num_cities(self) -> int:
        """Get the number of cities."""
        ...

    def edge_to_index(self, i: int, j: int) -> int:
        """
        Convert edge indices to variable index.

        Args:
            i: First city (i < j)
            j: Second city

        Returns:
            Variable index for edge (i,j)
        """
        ...

    def index_to_edge(self, idx: int) -> Tuple[int, int]:
        """
        Convert variable index to edge indices.

        Args:
            idx: Variable index

        Returns:
            Tuple (i, j) with i < j
        """
        ...
