add_executable(neutvect-converter neutvect-converter.cxx)

target_link_libraries(neutvect-converter PRIVATE nvconv)

install(TARGETS neutvect-converter EXPORT nvconv-targets)