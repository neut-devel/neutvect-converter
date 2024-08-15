cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(nvconv VERSION 0.9.0)

set(CMAKE_CXX_STANDARD 17)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}/"
           CACHE PATH "default install path" FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
elseif(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

find_package(ROOT 6 REQUIRED)
set(nvconv_MIN_NEUT_VERSION 5.5.0)
find_package(NEUT ${nvconv_MIN_NEUT_VERSION} REQUIRED)

find_package(HepMC3 3.3.0 REQUIRED)

include(CPM)
CPMAddPackage(
    NAME NuHepMC_CPPUtils
    GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
    GIT_TAG main
    OPTIONS "BUILTIN_HEPMC3 OFF"
)

add_subdirectory(src)
add_subdirectory(app)

install(EXPORT nvconv-targets
  FILE nvconvTargets.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/nvconv
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/nvconvConfig.cmake.in 
  ${CMAKE_BINARY_DIR}/nvconvConfig.cmake
  INSTALL_DESTINATION 
    /this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/nvconvConfigVersion.cmake
  VERSION ${nvconv_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES 
    ${CMAKE_BINARY_DIR}/nvconvConfig.cmake 
    ${CMAKE_BINARY_DIR}/nvconvConfigVersion.cmake 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/nvconv)