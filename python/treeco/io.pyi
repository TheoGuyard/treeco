from typing import List

from treeco.types import BinaryVector, RealVector
from treeco.geometry import Domain

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
