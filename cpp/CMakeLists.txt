cmake_minimum_required(VERSION 3.16)

project(treeco 
  VERSION 0.0.1 
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(TREECO_BUILD_EXEC "Build main executable" ON)
option(TREECO_BUILD_TESTS "Build tests" ON)
option(TREECO_BUILD_EXAMPLES "Build examples" ON)


# =========================================================================== #
# Find Gurobi library
# =========================================================================== #

if(NOT DEFINED ENV{GUROBI_HOME})
  message(FATAL_ERROR "Please ensure that Gurobi is installed and that the GUROBI_HOME environment variable is set.")
endif()

# Gurobi C++ header library
find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h
    HINTS ENV GUROBI_HOME
    PATH_SUFFIXES include
)

# Gurobi library (locate highest version number)
file(GLOB GUROBI_LIBRARY_CANDIDATES
  "$ENV{GUROBI_HOME}/lib/libgurobi*.so"
  "$ENV{GUROBI_HOME}/lib/libgurobi*.dylib"
)
set(GUROBI_LIBRARY "")
set(GUROBI_VERSION 0)
foreach(lib ${GUROBI_LIBRARY_CANDIDATES})
    string(REGEX MATCH "gurobi([0-9]+)" _ ${lib})
    if(CMAKE_MATCH_1)
        if(${CMAKE_MATCH_1} GREATER ${GUROBI_VERSION})
            set(GUROBI_VERSION ${CMAKE_MATCH_1})
            set(GUROBI_LIBRARY ${lib})
        endif()
    endif()
endforeach()
if(NOT GUROBI_LIBRARY)
    message(FATAL_ERROR "Failed to locate Gurobi library")
endif()

# Gurobi C++ library
find_library(GUROBI_CXX_LIBRARY
    NAMES gurobi_c++
    HINTS ENV GUROBI_HOME
    PATH_SUFFIXES lib
)

# Create imported targets for Gurobi library
add_library(gurobi::gurobi UNKNOWN IMPORTED GLOBAL)
set_target_properties(gurobi::gurobi PROPERTIES
    IMPORTED_LOCATION "${GUROBI_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
)

# Create imported targets for Gurobi C++ library
add_library(gurobi::gurobi_cxx UNKNOWN IMPORTED GLOBAL)
set_target_properties(gurobi::gurobi_cxx PROPERTIES
    IMPORTED_LOCATION "${GUROBI_CXX_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
)


# =========================================================================== #
# Treeco C++ library
# =========================================================================== #

add_library(treeco
  src/Adjacency.cpp
  src/CLI.cpp
  src/Dynprog.cpp
  src/Feasibility.cpp
  src/Geometry.cpp
  src/IO.cpp
  src/LDTree.cpp
  src/Problem.cpp
  src/Problem/Explicit.cpp
  src/Problem/Knapsack.cpp
  src/Problem/Tsp.cpp
  src/Problem/Maxcut.cpp
  src/Tree.cpp
  src/Voronoi.cpp
)

target_include_directories(treeco
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(treeco
  PUBLIC
    gurobi::gurobi
    gurobi::gurobi_cxx
)

target_compile_features(treeco PUBLIC cxx_std_20)

target_compile_options(treeco
  PRIVATE
    -Wall -Wextra -Wpedantic
)

# =========================================================================== #
# Treeco C++ executables
# =========================================================================== #

# Main executable
if(TREECO_BUILD_EXEC)
  add_executable(main src/main.cpp)
  target_link_libraries(main PRIVATE treeco)
endif()

# Tests executable
if(TREECO_BUILD_TESTS)
  enable_testing()

  include(FetchContent)
  
  message(STATUS "Fetching GTest...")
  FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG v1.14.0
  )
  FetchContent_MakeAvailable(googletest)

  add_executable(tests
    tests/test_adjacency.cpp
    tests/test_dynprog.cpp
    tests/test_feasibility.cpp
    tests/test_geometry.cpp
    tests/test_io.cpp
    tests/test_ldtree.cpp
    tests/test_problem.cpp
    tests/test_tree.cpp
    tests/test_voronoi.cpp
  )

  target_link_libraries(tests
    PRIVATE 
      treeco
      GTest::gtest
      GTest::gtest_main
  )
endif()

# Examples executable
if(TREECO_BUILD_EXAMPLES)
  add_executable(examples examples/example.cpp)
  target_link_libraries(examples PRIVATE treeco)
endif()
