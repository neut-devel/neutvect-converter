add_library(nvconv STATIC nvconv.cxx)

target_link_libraries(nvconv PUBLIC NEUT::IO HepMC3::All_static ROOT::RIO)
target_include_directories(nvconv PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_compile_definitions(nvconv PRIVATE 
  NEUT_VERSION_STR="5.5.0"
  PROJECT_VERSION_STR="${PROJECT_VERSION}")

set_target_properties(nvconv PROPERTIES 
  PUBLIC_HEADER nvconv.h)

install(TARGETS nvconv
    EXPORT nvconv-targets
    LIBRARY DESTINATION lib 
    ARCHIVE DESTINATION lib 
    PUBLIC_HEADER DESTINATION include)