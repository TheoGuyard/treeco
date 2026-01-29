"""Treeco - Linear Decision Tree Policies for Combinatorial Optimization."""

import sys

from . import _treeco

# Submodule imports
geometry = _treeco.geometry
problem = _treeco.problem
voronoi = _treeco.voronoi
dynprog = _treeco.dynprog
tree = _treeco.tree
ldtree = _treeco.ldtree
io = _treeco.io

# Submodule aliases
sys.modules[__name__ + ".geometry"] = geometry
sys.modules[__name__ + ".problem"] = problem
sys.modules[__name__ + ".voronoi"] = voronoi
sys.modules[__name__ + ".dynprog"] = dynprog
sys.modules[__name__ + ".tree"] = tree
sys.modules[__name__ + ".ldtree"] = ldtree
sys.modules[__name__ + ".io"] = io

# Top-level imports LDTree
LDTree = ldtree.LDTree

__all__ = [
    "LDTree",
    "geometry",
    "problem",
    "voronoi",
    "dynprog",
    "tree",
    "ldtree",
    "io",
]

__version__ = "0.0.1"
__authors__ = "Théo Guyard"
