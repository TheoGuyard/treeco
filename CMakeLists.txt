cmake_minimum_required(VERSION 3.16)
project(treeco VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Options
option(TREECO_BUILD_EXEC "Build executable" ON)
option(TREECO_BUILD_TESTS "Build tests" ON)
option(TREECO_BUILD_EXAMPLES "Build examples" ON)
option(TREECO_BUILD_PYTHON "Build Python bindings" OFF)


###### Third-party libraries

### Gurobi

if(NOT DEFINED ENV{GUROBI_HOME})
  message(FATAL_ERROR "Please set the GUROBI_HOME environment variable")
endif()

# Gurobi paths
set(GUROBI_HOME $ENV{GUROBI_HOME})
set(GUROBI_INCLUDE_DIRS "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIRS "${GUROBI_HOME}/lib")

# Gurobi library (locate highest version number)
file(GLOB GUROBI_CANDIDATES
  "${GUROBI_HOME}/lib/libgurobi*.so"
  "${GUROBI_HOME}/lib/libgurobi*.dylib"
)
set(GUROBI_LIBRARY "")
set(GUROBI_VERSION 0)
foreach(lib ${GUROBI_CANDIDATES})
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

# Gurobi C++ header library
find_library(
  GUROBI_CXX_LIBRARY
  NAMES gurobi_c++
  HINTS ${GUROBI_HOME}
  PATH_SUFFIXES lib
)


###### Main library

# Library sources
add_library(treeco_lib
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

# Library headers
target_include_directories(treeco_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link third-party libraries
target_include_directories(treeco_lib PUBLIC ${GUROBI_INCLUDE_DIRS})
target_link_directories(treeco_lib PUBLIC ${GUROBI_LIBRARY_DIRS})
target_link_libraries(treeco_lib PUBLIC
    ${GUROBI_LIBRARY}
    ${GUROBI_CXX_LIBRARY}
)


###### Main executable

if(TREECO_BUILD_EXEC)
  add_executable(treeco_exec src/main.cpp)
  target_link_libraries(treeco_exec PRIVATE treeco_lib)
  set_target_properties(treeco_exec PROPERTIES
      OUTPUT_NAME treeco
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()


###### Tests executable

if(TREECO_BUILD_TESTS)
    enable_testing()

    # Find testing framework
    find_package(GTest QUIET)

    # Fetch testing framework if not found
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
        )
        # For Windows: Prevent overriding the parent project's compiler/linker settings
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
        set(GTest_FOUND TRUE)
    endif()
      
    add_executable(treeco_tests
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
    
    target_link_libraries(treeco_tests PRIVATE
        treeco_lib
        GTest::gtest
        GTest::gtest_main
    )
    
    add_test(NAME treeco_tests COMMAND treeco_tests)

    set_target_properties(treeco_tests PROPERTIES
      OUTPUT_NAME tests
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()


###### Examples executable

if(TREECO_BUILD_EXAMPLES)
  add_executable(treeco_examples examples/example.cpp)
  target_link_libraries(treeco_examples PRIVATE treeco_lib)
  set_target_properties(treeco_examples PROPERTIES
      OUTPUT_NAME examples
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()


###### Python bindings

if(TREECO_BUILD_PYTHON)
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  
  # Find pybind11
  find_package(pybind11 CONFIG QUIET)
  
  # Fetch pybind11 if not found
  if(NOT pybind11_FOUND)
      include(FetchContent)
      FetchContent_Declare(
          pybind11
          GIT_REPOSITORY https://github.com/pybind/pybind11.git
          GIT_TAG v2.13.1
      )
      FetchContent_MakeAvailable(pybind11)
  endif()

  target_compile_definitions(treeco_lib PRIVATE TREECO_BUILD_PYTHON)
  
  pybind11_add_module(treeco_python MODULE
    python/treeco/_treeco.cpp
  )

  # match wrapper module name _treeco
  set_target_properties(treeco_python PROPERTIES
    OUTPUT_NAME _treeco
  )

  target_link_libraries(treeco_lib PRIVATE
    Python3::Python
    pybind11::headers
    pybind11::module
  )
  target_link_libraries(treeco_python PRIVATE treeco_lib)

  target_include_directories(treeco_python PRIVATE
      ${CMAKE_SOURCE_DIR}/include
      ${CMAKE_CURRENT_SOURCE_DIR}/treeco
  )
  install(TARGETS treeco_python LIBRARY DESTINATION treeco)
endif()


###### Installation

install(TARGETS treeco_lib
    EXPORT treecoTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/treeco DESTINATION include)

install(EXPORT treecoTargets
    FILE treecoTargets.cmake
    NAMESPACE treeco::
    DESTINATION lib/cmake/treeco
)
