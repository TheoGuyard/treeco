# Installation Guide

This document provides detailed instructions for building treeco from source and using its C++ library.

## Prerequisites

- **C++20** compatible compiler
- **CMake** 3.16 or higher
- **Gurobi** optimizer:
    1. Download and install Gurobi from [gurobi.com](https://www.gurobi.com/)
    2. Activate your license following Gurobi's instructions (free academic licenses are available)
    3. Set the `GUROBI_HOME` environment variable in your `~/.bashrc` or `~/.zshrc` file. It should look like:
       - For linux: `GUROBI_HOME=/opt/gurobi1100/linux64`
       - For macOS: `GUROBI_HOME=/Library/gurobi1100/macos_universal2`
    4. Restart your terminal


## Building the C++ Library

The C++ library can be downloaded and built using the following steps:

```bash
# Clone the repository
git clone https://github.com/TheoGuyard/treeco.git
cd treeco

# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make
```

Once done, the test and example executables can be run from the `build` directory as `./tests` and `./examples`.
The following options are available for the `cmake` command:

| Option | Default | Description |
|--------|---------|-------------|
| `TREECO_BUILD_EXEC`| `ON` | Build the main executable |
| `TREECO_BUILD_TESTS` | `ON` | Build tests executable |
| `TREECO_BUILD_EXAMPLES` | `ON` | Build examples executable |

and can be used as `cmake -DTREECO_BUILD_TESTS=OFF ..`.

## Using the Main Executable

The build process creates an executable in the `build` directory which can be run as follows:

```
./treeco --points <points_path> --domain <domain_path.txt> --queries <queries_path.txt> [options]
```

This executable builds a decision policy from a set of feasible points and a domain for the cost vector, then performs queries on the built tree.
The points, domain, and queries files must be formatted as follows:

```python
# points_path.txt
# - the first line give the dimension and number of points
# - the other lines represent each binary feasible points x in {0,1}^n
points n p
x11 x12 ... x1n
x21 x22 ... x2n
...
xp1 xp2 ... xpn

# domain_path.txt
# - the first line gives the dimension and number of inequalities modelling the cost domain
# - the other lines represent inequalities of the form ⟨a,x⟩ + b {LT,LE,EQ,GE,GT} 0 where LT = less than, LE = less equal, EQ = equal, GE = greater equal, GT = greater than.
domain n m
a11 a12 ... a1n b1 EQ
a21 a22 ... a2n b2 LE
...
am1 am2 ... amn bm GE

# queries_path.txt
# - the first line gives the dimension and number of queries
# - the other lines represent each real cost vector c in ℝ^n to be queried
queries n d
c11 c12 ... c1n
c21 c22 ... c2n
...
cd1 cd2 ... cdn
```

The domain file is optional and if not provided, the cost vectors are assumed to be unrestricted in ℝ^n.
Use `./treeco --help` to see all available options for the executable.

## Using the C++ Library

`treeco` C++ library can be used in your own projects by adding it.
Check out the [examples](examples) folder for a quick start guide.
The library can be liked to your `CMakeLists.txt` configuration file as follows:

```cmake
find_package(treeco REQUIRED)
target_link_libraries(your_target PRIVATE treeco::treeco)
```

Alternatively, you can directly include the `treeco` source code directory in your project and link against the library:

```cmake
add_subdirectory(treeco)
target_link_libraries(your_target PRIVATE treeco)
```
