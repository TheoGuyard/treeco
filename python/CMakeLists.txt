cmake_minimum_required(VERSION 3.16)

project(treeco 
  VERSION 0.0.1 
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("Fetching Python...")
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

message("Fetching Pybind11...")
find_package(pybind11 CONFIG REQUIRED)

# Force cpp disable unnecessary build options
set(TREECO_BUILD_EXEC OFF CACHE BOOL "" FORCE)
set(TREECO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(TREECO_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# Include cpp project
add_subdirectory(../cpp cpp_build)

# Link additional definition flags, headers, and libraries to treeco
# required when built for Python bindings (mostly to bind keyboard interrupts)
target_compile_definitions(treeco PRIVATE TREECO_BUILD_PYTHON)
target_include_directories(treeco PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(treeco PRIVATE pybind11::headers)

# Link pybind11 module
pybind11_add_module(_treeco treeco/_treeco.cpp)

# Link external libraries to the pybind11 module
target_link_libraries(_treeco PRIVATE treeco Python3::Python)

# Install pybind11 module
install(TARGETS _treeco LIBRARY DESTINATION treeco)