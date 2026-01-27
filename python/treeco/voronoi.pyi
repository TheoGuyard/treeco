from dataclasses import dataclass

@dataclass
class VoronoiStats:
    """
    Statistics from Voronoi diagram construction.

    Attributes:
        build_time: Total build time in seconds
        lp_solved: Number of LPs solved during construction
    """

    build_time: float
    lp_solved: int

class Voronoi:
    """
    Voronoi diagram for a feasible set of points.
    """

    @property
    def dim_points(self) -> int:
        """Get the dimension of the points."""
        ...

    @property
    def num_points(self) -> int:
        """Get the number of feasible points."""
        ...

    @property
    def num_splits(self) -> int:
        """Get the number of Voronoi splits."""
        ...

    @property
    def num_faces(self) -> int:
        """Get the number of Voronoi faces."""
        ...

    @property
    def num_edges(self) -> int:
        """Get the number of Voronoi edges."""
        ...

    @property
    def stats(self) -> VoronoiStats:
        """Get construction statistics."""
        ...
