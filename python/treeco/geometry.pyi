from enum import Enum
from typing import List, Tuple

from treeco.types import RealVector

class Relation(Enum):
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

    lt: Relation
    le: Relation
    eq: Relation
    ge: Relation
    gt: Relation
    rt: Relation
    rf: Relation

ConstrData = Tuple[RealVector, float, Relation]
"""Tuple (a,b,rel) for linear constraint <a,x> + b rel 0."""

Domain = List[ConstrData]
"""Collection of linear constraints defining a polyhedral domain."""

@staticmethod
def positive_orthant(dim: int) -> Domain:
    """
    Generate the positive orthant domain in R^dim.

    Args:
        dim: Dimension of the space
    Returns:
        Domain representing x_i > 0 for all i
    """
    ...

@staticmethod
def non_negative_orthant(dim: int) -> Domain:
    """
    Generate the non-negative orthant domain in R^dim.

    Args:
        dim: Dimension of the space
    Returns:
        Domain representing x_i >= 0 for all i
    """
    ...

@staticmethod
def negative_orthant(dim: int) -> Domain:
    """
    Generate the negative orthant domain in R^dim.

    Args:
        dim: Dimension of the space
    Returns:
        Domain representing x_i < 0 for all i
    """
    ...

@staticmethod
def non_positive_orthant(dim: int) -> Domain:
    """
    Generate the non-positive orthant domain in R^dim.

    Args:
        dim: Dimension of the space
    Returns:
        Domain representing x_i <= 0 for all i
    """
    ...
