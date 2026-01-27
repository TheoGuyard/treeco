# Treeco - C++ Interface

[![C++20](https://img.shields.io/badge/C++-20-blue.svg)](https://isocpp.org/)

The C++ interface to `treeco` can be downloaded as follows:

```bash
# Clone the repository
git clone https://github.com/TheoGuyard/treeco.git
```

It can be built from source from the `cpp/` directory as follows:

```bash
cd treeco/cpp               # Navigate to the cpp folder
mkdir build && cd build     # Create build directory
cmake ..                    # Configure
make                        # Build
```

The following options are available for the `cmake` command during build:

| Option                  | Default | Description               |
| ----------------------- | ------- | ------------------------- |
| `TREECO_BUILD_EXEC`     | `ON`    | Build the main executable |
| `TREECO_BUILD_TESTS`    | `ON`    | Build tests executable    |
| `TREECO_BUILD_EXAMPLES` | `ON`    | Build examples executable |

and can be used as `cmake -DTREECO_BUILD_TESTS=OFF ..`.

After building, `treeco` library is available in the `cpp/lib` folder. All executables are located in the `cpp/bin` folder:

* `./treeco` to call the library from the command line (use the `--help` flag to obtain usage instructions).
* `./tests` to run unit tests
* `./examples` to run the example code

Check out the [cpp/examples](cpp/examples) folder for a quick start guide. 

> A more detailed documentation is coming soon!