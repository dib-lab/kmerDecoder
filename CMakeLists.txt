# ----------------------------------------------------------------------------
# kmerDecoder project
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Project Info ----

project(
        kmerDecoder
        VERSION 1.0
        LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(
            FATAL_ERROR
            "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
    )
endif()

# ---- Add dependencies via CPM ----

include(cmake/CPM.cmake)

# PackageProject.cmake will be used to make our target installable
CPMAddPackage(
        NAME PackageProject.cmake
        GITHUB_REPOSITORY TheLartians/PackageProject.cmake
        VERSION 1.4.1
)


# ---- Create binary ----

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_BIN "Build Binaries" OFF)

if (BUILD_BIN)
    add_definitions(-BUILD_BIN=1)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fPIC -fopenmp -pedantic -O3 -oFast")

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/src/kmerDecoder.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_kmers.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_skipmers.cpp"
        "${PROJECT_SOURCE_DIR}/src/KD_minimizers.cpp"
        "${PROJECT_SOURCE_DIR}/src/HashUtils/hashutil.cpp"
        "${PROJECT_SOURCE_DIR}/src/HashUtils/aaHasher.cpp"
        "${PROJECT_SOURCE_DIR}/src/Utils/kmer.cpp"
        )

# TO DO delete when no need
# file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*h" "${CMAKE_CURRENT_SOURCE_DIR}/include/*hpp")
# file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")


set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap")
set(DOCTEST_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/doctest/")
set(KMERDECODER_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(KSEQ_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/kseq")

# include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(${PHMAP_INCLUDE_DIRS})
# include_directories(${DOCTEST_INCLUDE_DIRS})
include_directories(${KSEQ_INCLUDE_DIRS})

# ---- Create library ----

add_library(kmerDecoder ${SOURCES} ${KMERDECODER_INCLUDE_DIRS} ${KSEQ_INCLUDE_DIRS} ${PHMAP_INCLUDE_DIRS})

set_target_properties(kmerDecoder PROPERTIES CXX_STANDARD 17)

target_include_directories(
        kmerDecoder PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

#target_include_directories(kmerDecoder INTERFACE ${PHMAP_INCLUDE_DIRS} ${KSEQ_INCLUDE_DIRS} ${KMERDECODER_INCLUDE_DIRS})
# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(kmerDecoder PUBLIC "$<$<BOOL:${MSVC}>:/permissive->")
target_link_libraries(kmerDecoder z)

if (BUILD_BIN)
    message("Building kmerDecoder binaries")
    add_subdirectory(bin)
endif ()


# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)


# ---- Package Info -------

packageProject(
        NAME ${PROJECT_NAME}
        VERSION ${PROJECT_VERSION}
        NAMESPACE ${PROJECT_NAME}
        BINARY_DIR ${PROJECT_BINARY_DIR}
        INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
        INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
        VERSION_HEADER "${VERSION_HEADER_LOCATION}"
        DEPENDENCIES ""
)