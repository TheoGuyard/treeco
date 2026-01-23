# Treeco

[![License: MIT](https://img.shields.io/badge/License-AGPL--v3-red.svg)](https://github.com/TheoGuyard/treeco/blob/main/LICENSE)
[![C++20](https://img.shields.io/badge/C++-20-blue.svg)](https://isocpp.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)


`treeco` builds efficient decision policies for combinatorial optimization problems of the form:

```
max  ⟨c, x⟩   s.t.  x ∈ X
```

where `X ⊆ {0,1}^n` is a feasible set and `c ∈ ℝ^n` is a cost vector. The policy is **built once** to encode a **fixed** feasible set into a linear decision tree, then can be **queried** efficiently to retrieve optimal solutions for **any** cost vector by traversing the tree.

## Installation

`treeco` can be installed using [pip](https://pypi.org/project/pip/) with the command
```bash
pip install git+https://github.com/TheoGuyard/treeco.git#subdirectory=python
```
Check out the [python/examples](python/examples) folder for a quick start guide. A more detailed documentation is coming soon!

**Requirements:** The installation requires [Gurobi](https://www.gurobi.com/) solver with the environment variable `GUROBI_HOME` properly set to its main directory path.
Future versions will support open-source solvers to remove this dependency.

> For C++ library installation, tests, examples, and advanced options, see [INSTALL.md](INSTALL.md).

## Troubleshooting

For any troubleshooting during installation, please open an [issue](https://github.com/TheoGuyard/treeco/issues) on the Github repository.

## License

Treeco is distributed under the [AGPL v3 license](LICENSE).
