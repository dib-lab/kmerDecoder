# Minimum cmake version
cmake_minimum_required (VERSION 3.4)

project (kmerDecoder C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")


file(GLOB SOURCES 
"src/kmers.cpp" 
"src/skipmers.cpp"
"src/minimizers.cpp"
)

include_directories("include")

add_library(kmerDecoder STATIC ${SOURCES})