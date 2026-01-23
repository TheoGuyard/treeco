"""
treeco - Linear Decision Tree Policies for Combinatorial Optimization

This module provides Python bindings for the treeco C++ library, enabling
efficient decision policies for combinatorial optimization problems.
"""

from typing import List, Any, Tuple, overload

# Type aliases
BinaryVector = List[int]
"""Binary vector in {0,1}^n representing a feasible solution."""

RealVector = List[float]
"""Real-valued vector in R^n representing a cost vector."""

ConstrData = Tuple[RealVector, float, "Relation"]
"""Tuple (a,b,rel) for linear constraint <a,x> + b rel 0."""

Domain = List[ConstrData]
"""Collection of linear constraints defining a polyhedral domain."""


class Exploration:
    """
    Exploration strategy for the dynamic programming algorithm.
    
    Attributes:
        greedy: Single iteration with k=1 (fast, potentially suboptimal)
        iterative: Iterative scheme with k=1,2,...,numSplits (optimal, early stopping)
        exhaustive: Single iteration with k=numSplits (optimal, no early stopping)
    """
    greedy: "Exploration"
    iterative: "Exploration"
    exhaustive: "Exploration"


class Branching:
    """
    Branching mode for decision tree nodes.
    
    Attributes:
        ternary: Three-way branching: {<, =, >}
        binary: Two-way branching: {<, >} with tie-breaking for {=}
    """
    ternary: "Branching"
    binary: "Branching"


class LowerBounding:
    """
    Lower bound computation strategy for subtree depth.
    
    Attributes:
        fixed: Fixed lower bound: ceil(log_c(|F|)) where c is child count
        backtrack: Backtrack from children's bounds for tighter estimate
    """
    fixed: "LowerBounding"
    backtrack: "LowerBounding"


class Positioning:
    """
    Position computation mode for face-split relations.
    
    Attributes:
        online: Compute positions on-demand during search
        precompute: Precompute all face-split positions upfront
    """
    online: "Positioning"
    precompute: "Positioning"


class SplitSelection:
    """
    Split selection strategy at each decision node.
    
    Attributes:
        all: Consider all valid splits at each node
        sampling: Randomly sample k splits to score
    """
    all: "SplitSelection"
    sampling: "SplitSelection"


class SplitScoring:
    """
    Split scoring strategy for ordering candidate splits.
    
    Attributes:
        variance: Variance-based deviation from equal split
        entropy: Information gain (entropy reduction)
        minmax: Minimize maximum child face count
        none: All splits have equal score
        random: Random scoring
    """
    variance: "SplitScoring"
    entropy: "SplitScoring"
    minmax: "SplitScoring"
    none: "SplitScoring"
    random: "SplitScoring"


class Relation:
    """
    Relation type for linear constraints.
    
    Represents the position of a point relative to a hyperplane <a,x> + b.
    
    Attributes:
        lt: Strict less than (<)
        le: Less than or equal (<=)
        eq: Equal (=)
        ge: Greater than or equal (>=)
        gt: Strict greater than (>)
        rt: Always true (unconstrained)
        rf: Always false (infeasible)
    """
    lt: "Relation"
    le: "Relation"
    eq: "Relation"
    ge: "Relation"
    gt: "Relation"
    rt: "Relation"
    rf: "Relation"


class LDTreeStats:
    """
    Statistics from LDTree construction.
    
    Attributes:
        build_time: Total build time in seconds
    """
    build_time: float
    
    def __init__(self) -> None:
        """Initialize empty statistics."""
        ...


class LDTree:
    """
    Linear Decision Tree for combinatorial optimization.
    
    This class provides a complete workflow for building decision trees that
    can efficiently solve linear optimization problems over a fixed feasible set.
    Once built, the tree can be queried with any cost vector to retrieve the
    optimal solutions in O(depth) time.
    
    Example:
        >>> tree = LDTree(points)
        >>> tree.build()
        >>> solutions = tree.query(cost)
    """
    
    @overload
    def __init__(self, file_points: str, file_domain: str = "") -> None:
        """
        Construct an LDTree from input files.
        
        Args:
            file_points: Path to file containing feasible points
            file_domain: Path to file containing cost domain constraints (optional)
        """
        ...
    
    @overload
    def __init__(self, points: List[BinaryVector], domain: Domain = ...) -> None:
        """
        Construct an LDTree from data objects.
        
        Args:
            points: List of feasible binary solutions
            domain: Optional linear constraints on the cost domain
        """
        ...
    
    def build(
        self,
        verbose: bool = True,
        log_interval: float = 5.0,
        time_limit: float = ...,
        tolerance: float = 1e-8,
        deduplicate: bool = True,
        filter_checks: bool = True,
        exploration: Exploration = ...,
        branching: Branching = ...,
        lower_bounding: LowerBounding = ...,
        positioning: Positioning = ...,
        split_selection: SplitSelection = ...,
        split_scoring: SplitScoring = ...,
        random_seed: int = 42
    ) -> None:
        """
        Build the LDTree structure.
        
        Constructs the Voronoi diagram and decision tree structure using
        dynamic programming to find the minimum-depth tree.
        
        Args:
            verbose: Enable verbose logging output
            log_interval: Logging interval in seconds
            time_limit: Maximum build time in seconds
            tolerance: Numerical tolerance for equality comparisons
            deduplicate: Remove duplicate feasible points before building
            filter_checks: Use filtering for faster split validity checks
            exploration: Search exploration strategy
            branching: Tree branching mode (binary or ternary)
            lower_bounding: Lower bound computation strategy
            positioning: Face-split position computation mode
            split_selection: Split candidate selection strategy
            split_scoring: Split quality scoring strategy
            random_seed: Random seed for sampling-based selection
        """
        ...
    
    def query(self, cost: RealVector, check_domain: bool = False) -> List[BinaryVector]:
        """
        Query the tree for optimal solutions.
        
        Args:
            cost: The cost vector to optimize
            check_domain: Validate that cost is within the domain before querying
            
        Returns:
            List of optimal binary solutions
        """
        ...
    
    def pprint(self, tight_display: bool = False) -> None:
        """
        Pretty-print the tree structure.
        
        Args:
            tight_display: If True, show only indices; otherwise show full vectors
        """
        ...
    
    def flatten(self, filepath: str, doc: str = "", benchmark_mode: bool = False) -> None:
        """
        Generate standalone C code implementing the decision tree.
        
        Args:
            filepath: Output file path for the generated code
            doc: Documentation header for the generated code
            benchmark_mode: Generate benchmarking code instead of query code
        """
        ...
    
    def domain(self) -> Domain:
        """Get the cost domain constraints."""
        ...
    
    def voronoi(self) -> Any:
        """Get the Voronoi diagram."""
        ...
    
    def tree(self) -> Any:
        """Get the decision tree structure."""
        ...
    
    def stats(self) -> LDTreeStats:
        """Get build statistics."""
        ...


class Problem:
    """
    Abstract base class for optimization problems.
    
    Represents a linear optimization problem of the form:
        max <c, x>  subject to  x in X
    
    where c is a real-valued cost vector and X is a set of binary vectors.
    """
    
    def dimension(self) -> int:
        """Get the dimension of the problem (length of vectors x and c)."""
        ...
    
    def get_feasible_set(self) -> List[BinaryVector]:
        """Get all feasible solutions as a list of binary vectors."""
        ...
    
    def sample_cost(self) -> RealVector:
        """Sample a random cost vector with entries uniformly in [0, 1]."""
        ...


class Explicit(Problem):
    """
    Optimization problem with an explicitly provided feasible set.
    
    Use this class when you have a pre-computed list of all feasible solutions.
    """
    
    def __init__(self, feasible_set: List[BinaryVector]) -> None:
        """
        Construct from an explicit feasible set.
        
        Args:
            feasible_set: List of all feasible binary solutions
        """
        ...
    
    def get_feasible_set(self) -> List[BinaryVector]:
        """Get all feasible solutions."""
        ...


class Maxcut(Problem):
    """
    Maximum Cut Problem on a complete undirected graph.
    
    Given d vertices, partition vertices into two sets S and T such that
    the total weight of edges crossing the partition is maximized.
    """
    
    def __init__(self, num_vertices: int) -> None:
        """
        Construct a MaxCut problem.
        
        Args:
            num_vertices: Number of vertices in the complete graph
        """
        ...
    
    def num_vertices(self) -> int:
        """Get the number of vertices."""
        ...
    
    def get_feasible_set(self) -> List[BinaryVector]:
        """Get all feasible solutions (all valid cuts)."""
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
    
    Given d cities, find a Hamiltonian cycle (tour visiting each city
    exactly once) that minimizes the total distance.
    """
    
    def __init__(self, num_cities: int) -> None:
        """
        Construct a TSP instance.
        
        Args:
            num_cities: Number of cities in the complete graph
        """
        ...
    
    def num_cities(self) -> int:
        """Get the number of cities."""
        ...
    
    def get_feasible_set(self) -> List[BinaryVector]:
        """Get all feasible solutions (all valid Hamiltonian cycles)."""
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


# IO functions

def read_points(filepath: str) -> List[BinaryVector]:
    """
    Read feasible points from a file.
    
    Args:
        filepath: Path to the points file
        
    Returns:
        List of binary vectors representing feasible solutions
    """
    ...


def write_points(filepath: str, points: List[BinaryVector]) -> None:
    """
    Write feasible points to a file.
    
    Args:
        filepath: Output file path
        points: List of binary vectors to write
    """
    ...


def read_domain(filepath: str) -> Domain:
    """
    Read domain constraints from a file.
    
    Args:
        filepath: Path to the domain file
        
    Returns:
        List of constraint tuples (coefficients, relation, rhs)
    """
    ...


def write_domain(filepath: str, domain: Domain) -> None:
    """
    Write domain constraints to a file.
    
    Args:
        filepath: Output file path
        domain: List of constraint tuples to write
    """
    ...


def read_queries(filepath: str) -> List[RealVector]:
    """
    Read query cost vectors from a file.
    
    Args:
        filepath: Path to the queries file
        
    Returns:
        List of real-valued cost vectors
    """
    ...


def write_queries(filepath: str, queries: List[RealVector]) -> None:
    """
    Write query cost vectors to a file.
    
    Args:
        filepath: Output file path
        queries: List of cost vectors to write
    """
    ...
