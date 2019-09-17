# ----------------------------------------------------------------------------
# kmerDecoder project
# ----------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.4)
project (kmerDecoder C CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")



# ----------------------------------------------------------------------------
# Seqan Setup
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
set(CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/lib/seqan/util/cmake")
set(SEQAN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/seqan/include")
message("${PROJECT_SOURCE_DIR}/lib/seqan/include")
find_package (SeqAn REQUIRED)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
message("${SEQAN_INCLUDE_DIRS}")
# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# ----------------------------------------------------------------------------
# END Seqan Setup
# ----------------------------------------------------------------------------

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -fPIC -fopenmp -W -pedantic -O3 -oFast")


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/src/kmerDecoder.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_kmers.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_skipmers.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_minimizers.cpp"
        "${PROJECT_SOURCE_DIR}/src/HashUtils/hashutil.cpp"
        "${PROJECT_SOURCE_DIR}/src/Utils/kmer.cpp"
)

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap")
set(KMERDECODER_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(${PHMAP_INCLUDE_DIRS})


add_library(kmerDecoder STATIC ${SOURCES})
target_include_directories(kmerDecoder INTERFACE ${PHMAP_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/lib/seqan/include ${KMERDECODER_INCLUDE_DIRS})
target_link_libraries(kmerDecoder z)
