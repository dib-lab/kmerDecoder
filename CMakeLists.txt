# ----------------------------------------------------------------------------
# kmerDecoder project
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

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
        VERSION 1.6.0
)


# ---- Create binary ----

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_BIN "Build Binaries" OFF)

if (BUILD_BIN)
    add_definitions(-BUILD_BIN=1)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fPIC -fopenmp -O3 -Ofast")


# ---- Create library ----

CPMAddPackage(
        NAME kseq
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/lib/kseq
)

CPMAddPackage(
        NAME PHMAP
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/lib/parallel-hashmap
)

if (PHMAP_ADDED)

 add_library(PHMAP INTERFACE IMPORTED)
 target_include_directories(PHMAP INTERFACE ${PHMAP_SOURCE_DIR})

endif()


set(SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/KD_kmers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/kmerDecoder.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/KD_skipmers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/KD_minimizers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/KD_aakmers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/hashUtils/aaHasher.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/hashUtils/hashutil.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/kmer.cpp"
        )

add_library(kmerDecoder ${SOURCES})

set_target_properties(kmerDecoder PROPERTIES CXX_STANDARD 14)

set(PHMAP_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/lib/parallel-hashmap")
set(KMERDECODER_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
set(KSEQ_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/lib/kseq")


target_include_directories(
        kmerDecoder PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

#target_include_directories(kmerDecoder INTERFACE ${PHMAP_INCLUDE_DIRS} ${KSEQ_INCLUDE_DIRS} ${KMERDECODER_INCLUDE_DIRS})
# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(kmerDecoder PUBLIC "$<$<BOOL:${MSVC}>:/permissive->")
target_link_libraries(kmerDecoder z kseq::kseq PHMAP)

if (BUILD_BIN)
    message("Building kmerDecoder binaries")
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
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
        DEPENDENCIES "kseq 1.3; PHMAP 1.0.0"
)