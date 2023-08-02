cmake_minimum_required(VERSION 3.22)
project(DeGenPrime)
add_executable(degenprime datanode.cpp  datasequence.cpp  degenprime.cpp  format.cpp  globalsettings.cpp  primercalculator.cpp  primer.cpp  primerpair.cpp  primerpairlist.cpp  sequence.cpp  sequencelist.cpp  sequencereader.cpp)
install(TARGETS degenprime DESTINATION bin)