add_library(nvconv STATIC nvconv.cxx)

target_link_libraries(nvconv PUBLIC NEUT::IO HepMC3::All ROOT::RIO)
target_include_directories(nvconv PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_target_properties(nvconv PROPERTIES 
  PUBLIC_HEADER nvconv.h)

install(TARGETS nvconv
    EXPORT nvconv-targets
    LIBRARY DESTINATION lib 
    ARCHIVE DESTINATION lib 
    PUBLIC_HEADER DESTINATION include)