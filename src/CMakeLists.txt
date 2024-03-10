cmake_minimum_required(VERSION 3.11)

# This is the version ID of degenprime.
project(DeGenPrime VERSION 0.1.2 LANGUAGES CXX)

configure_file(config.h.in config.h)

add_executable(degenprime datanode.cpp datasequence.cpp degenprime.cpp format.cpp globalsettings.cpp primercalculator.cpp primer.cpp primerpair.cpp primerpairlist.cpp sequence.cpp sequencelist.cpp sequencereader.cpp)
target_compile_options(degenprime PUBLIC -Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused)
set_property(TARGET degenprime PROPERTY CXX_STANDARD 17)
set_property(TARGET degenprime PROPERTY CXX_STANDARD_REQUIRED true)
set_property(TARGET degenprime PROPERTY DEBUG_POSTFIX -debug)

target_include_directories(degenprime PUBLIC "${PROJECT_BINARY_DIR}")

install(TARGETS degenprime DESTINATION bin)

set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
