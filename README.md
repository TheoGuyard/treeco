# Treeco

`treeco` is a library to build efficient decision policies for combinatorial optimization problems of the form:

```
max ⟨c,x⟩ s.t. x ∈ X
```

where `X ⊆ {0,1}^n` is a binary feasible set and `c ∈ C` is a cost vector within some domain `C ⊆ ℝ^n`. The policy is **built once** to encode a **fixed** feasible set into a linear decision tree, and can then be **queried** efficiently to retrieve optimal solutions corresponding to **any** cost vector by traversing the tree.

Check out the [C++ interface](./cpp/README.md) or the [Python interface](./python/README.md) for more details on how to install and use the library.


> **Requirements:** `treeco` currently requires [Gurobi](https://www.gurobi.com/) solver to be installed, with the environment variable `GUROBI_HOME` properly set to its main directory path. Future versions will support open-source solvers to remove this dependency.


## IO Format

Besides direct access to `treeco` internal functionalities to define and manipulate optimization problems, data for the feasible set, the cost vector domain, and the cost vector queries can be provided to the `treeco` library via `.txt` files in a specific format.

#### Feasible set
Gives a set of binary feasible points `X ⊆ {0,1}^n`. The header line contains the  keyword `points`, the dimension, and number of points. Each subsequent line contains the data of one feasible point. An example of a feasible set file with 3 points in `{0,1}^4` is as follows:

```python
# points.txt
points 4 3
0 0 1 1  # x1
0 1 0 1  # x2
1 0 0 0  # x3
```

#### Cost vector domain
Gives the cost vector domain `C ⊆ R^n` modelled as a set of linear inequalities. The header line contains the keyword `domain`, the dimension, and number of inequalities. Each subsequent line contains the coefficients of one inequality in the form `⟨a,x⟩ + b {LT,LE,EQ,GE,GT} 0` where `LT = less than, LE = less equal, EQ = equal, GE = greater equal, GT = greater than`. An example of a cost vector domain file in `ℝ^3` is as follows:

```python
# domain.txt
domain 4 3
 1.0  2.0 -1.0  0.0  0.5  EQ  #  x1 + 2x2 - x3 + 0.5  = 0
-1.0  0.0  1.0  0.0  1.0  GE  # -x1       + x3 + 1   >= 0
 0.0  1.0  1.0  0.0 -2.0  LE  #        x2 + x3 - 2   <= 0
```

#### Cost vector queries
Gives a pool of cost vectors to be queried. The header line contains the keyword `queries`, the dimension, and number of queries. Each subsequent line contains one cost vector. An example of a cost vector queries file with 3 queries in `ℝ^4` is as follows:

```python
# queries.txt
queries 4 3
 1.0  0.5 -1.0  2.0  # c1
-0.5  2.0  0.0  1.0  # c2
 0.0 -1.0  1.0  0.5  # c3
```

## Troubleshooting

For any troubleshooting during installation, please open an [issue](https://github.com/TheoGuyard/treeco/issues) on the Github repository.

## Licence

`treeco` is released under the [AGPL-3.0](LICENSE) licence.