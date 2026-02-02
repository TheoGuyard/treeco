from dataclasses import dataclass
from enum import Enum
from typing import List

class Exploration(Enum):
    """
    Exploration strategy for the dynamic programming algorithm.

    Attributes:
        greedy: Single iteration with k=1 (fast, potentially suboptimal)
        iterative: Iterative scheme with k=1,2,...,numSplits (optimal, early stopping)
        exhaustive: Single iteration with k=numSplits (optimal, no early stopping)
    """

    greedy: Exploration
    iterative: Exploration
    exhaustive: Exploration

class Branching(Enum):
    """
    Branching mode for decision tree nodes.

    Attributes:
        ternary: Three-way branching: {<, =, >}
        binary: Two-way branching: {<, >} with tie-breaking for {=}
    """

    ternary: Branching
    binary: Branching

class SplitScoring(Enum):
    """
    Split scoring strategy for ordering candidate splits.

    Attributes:
        variance: Variance-based deviation from equal split
        entropy: Information gain (entropy reduction)
        minmax: Minimize maximum child face count
        none: All splits have equal score
        random: Random scoring
    """

    variance: SplitScoring
    entropy: SplitScoring
    minmax: SplitScoring
    none: SplitScoring
    random: SplitScoring

@dataclass
class DynprogStats:
    """
    Statistics from dynamic programming tree construction.

    Attributes:
        run_time: Total run time in seconds
        num_iters: Number of DP iterations performed
        num_evals: Number of split evaluations
        num_states: Total number of DP states processed
        num_states_built: Number of DP states built
        num_states_closed: Number of DP states closed
        num_states_leafed: Number of DP states leafed
        num_states_pruned: Number of DP states pruned
        lp_solved: Number of feasibility checks performed
        optimal_depth: Optimal tree depth found
    """

    run_time: float
    num_iters: int
    num_evals: int
    num_states: int
    num_states_built: int
    num_states_closed: int
    num_states_leafed: int
    num_states_pruned: int
    lp_solved: int
    optimal_depth: int

DynprogLogs = List[DynprogStats]
