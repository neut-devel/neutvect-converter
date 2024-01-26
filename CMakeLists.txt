cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(neutvect-converter VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT 6 REQUIRED)
find_package(NEUT 5.5.0 REQUIRED)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

include(CPM)
CPMFindPackage(
    NAME NuHepMC_CPPUtils
    GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
    GIT_TAG main
    OPTIONS "BUILTIN_HEPMC3 ON"
            "NuHepMC_BUILD_STATIC_LIBS ON"
)

add_subdirectory(src)
add_subdirectory(app)