from dataclasses import dataclass
from typing import List, overload

from treeco.types import BinaryVector, RealVector
from treeco.geometry import Domain
from treeco.voronoi import Voronoi
from treeco.dynprog import (
    Exploration,
    Branching,
    LowerBounding,
    Positioning,
    SplitSelection,
    SplitScoring,
)
from treeco.tree import Tree

@dataclass
class LDTreeStats:
    """
    Statistics from LDTree construction.

    Attributes:
        build_time: Total build time in seconds
    """

    build_time: float

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
    def __init__(
        self,
        points: List[BinaryVector],
        domain: Domain = ...,
    ) -> None:
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
        log_save: bool = True,
        time_limit: float = float("inf"),
        tolerance: float = 1e-8,
        deduplicate: bool = True,
        filter_checks: bool = True,
        exploration: Exploration = Exploration.iterative,
        branching: Branching = Branching.binary,
        lower_bounding: LowerBounding = LowerBounding.backtrack,
        positioning: Positioning = Positioning.precompute,
        split_selection: SplitSelection = SplitSelection.all,
        split_scoring: SplitScoring = SplitScoring.variance,
        random_seed: int = 42,
    ) -> None:
        """
        Build the LDTree structure.

        Constructs the Voronoi diagram and decision tree structure using
        dynamic programming to find the minimum-depth tree.

        Args:
            verbose: Enable verbose logging output
            log_interval: Logging interval in seconds
            log_save: Save logs during dynamic programming
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

    def query(
        self, cost: RealVector, check_domain: bool = False
    ) -> List[BinaryVector]:
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

    def flatten(
        self, filepath: str, doc: str = "", benchmark_mode: bool = False
    ) -> None:
        """
        Generate standalone C code implementing the decision tree.

        Args:
            filepath: Output file path for the generated code
            doc: Documentation header for the generated code
            benchmark_mode: Generate benchmarking code instead of query code
        """
        ...

    @property
    def domain(self) -> Domain:
        """Get the cost domain constraints."""
        ...

    @property
    def voronoi(self) -> Voronoi:
        """Get the Voronoi diagram."""
        ...

    @property
    def tree(self) -> Tree:
        """Get the decision tree structure."""
        ...

    @property
    def stats(self) -> LDTreeStats:
        """Get build statistics."""
        ...
