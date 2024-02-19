add_library(nvconv SHARED nvconv.cxx nvfatxtools.cxx)

target_link_libraries(nvconv PUBLIC NEUT::IO NuHepMC::CPPUtils ROOT::RIO)
target_include_directories(nvconv PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(nvconv PRIVATE NEUTCONV_DEBUG)
endif()

target_compile_definitions(nvconv PRIVATE 
  NEUT_VERSION_STR="${NEUT_VERSION}"
  PROJECT_VERSION_STR="${PROJECT_VERSION}")

set_target_properties(nvconv PROPERTIES 
  PUBLIC_HEADER "nvconv.h;nvfatxtools.h")

install(TARGETS nvconv
    EXPORT nvconv-targets
    LIBRARY DESTINATION lib 
    ARCHIVE DESTINATION lib 
    PUBLIC_HEADER DESTINATION include)