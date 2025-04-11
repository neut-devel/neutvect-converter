add_executable(neutvect-converter neutvect-converter.cxx)

target_link_libraries(neutvect-converter PRIVATE nvconv)

set_target_properties(neutvect-converter PROPERTIES 
  INSTALL_RPATH "\${ORIGIN}/../lib")

install(TARGETS neutvect-converter EXPORT nvconv-targets)