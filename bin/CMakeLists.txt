cmake_minimum_required (VERSION 3.4)

project (kmerDecoder C CXX)

add_executable(test_kmerDecoder main.cpp)

file(COPY sample.fa DESTINATION ${CMAKE_BINARY_DIR})
file(COPY aaSample.fa DESTINATION ${CMAKE_BINARY_DIR})


target_link_libraries(test_kmerDecoder kmerDecoder)