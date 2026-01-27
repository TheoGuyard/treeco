from dataclasses import dataclass

from treeco.dynprog import DynprogStats

@dataclass
class TreeStats:
    """
    Statistics from LDTree construction.

    Attributes:
        build_time: Total build time in seconds
    """

    build_time: float
    dynprog_stats: DynprogStats

class Tree:
    """
    Tree structure
    """

    @property
    def size(self) -> int:
        """Get the size of the tree (number of nodes)."""
        ...

    @property
    def width(self) -> int:
        """Get the width of the tree."""
        ...

    @property
    def depth(self) -> int:
        """Get the depth of the tree."""
        ...

    @property
    def stats(self) -> TreeStats:
        """Get build statistics."""
        ...
