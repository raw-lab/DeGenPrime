cmake_minimum_required(VERSION 3.22)
project(DeGenPrime)
add_executable(DeGenPrime DataNode.cpp  DataSequence.cpp  DeGenPrime.cpp  GlobalSettings.cpp  PrimerCalculator.cpp  Primer.cpp  PrimerPair.cpp  PrimerPairList.cpp  Sequence.cpp  SequenceList.cpp  SequenceReader.cpp)
install(TARGETS DeGenPrime DESTINATION bin)