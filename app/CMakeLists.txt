add_executable(neutvect-converter neutvect-converter.cxx)

target_link_libraries(neutvect-converter PRIVATE nvconv)

add_executable(DumpNuHepMCSplines DumpNuHepMCSplines.cxx)
target_link_libraries(DumpNuHepMCSplines PRIVATE NuHepMC::CPPUtils_static ROOT::Hist ROOT::RIO)

install(TARGETS neutvect-converter DumpNuHepMCSplines EXPORT nvconv-targets)