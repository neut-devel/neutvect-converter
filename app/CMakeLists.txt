add_executable(neutvect-converter neutvect-converter.cxx)

target_link_libraries(neutvect-converter PRIVATE nvconv)

if(ZLIB_FOUND)
  target_link_libraries(neutvect-converter PRIVATE ZLIB::ZLIB)
  target_compile_definitions(neutvect-converter PRIVATE 
    HEPMC3_USE_COMPRESSION
    HEPMC3_Z_SUPPORT)
endif()


install(TARGETS neutvect-converter EXPORT nvconv-targets)