cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(kmerDecoderTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

CPMAddPackage(
        NAME doctest
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/../lib/doctest
)

if(TEST_INSTALLED_VERSION)
    find_package(kmerDecoder REQUIRED)
else()
    CPMAddPackage(NAME kmerDecoder SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

CPMAddPackage(
        NAME Format.cmake
        GITHUB_REPOSITORY TheLartians/Format.cmake
        VERSION 1.6
        OPTIONS # enable cmake formatting
        "FORMAT_CHECK_CMAKE ON"
)

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp)
add_executable(kmerDecoderTests ${sources})
target_link_libraries(kmerDecoderTests doctest kmerDecoder::kmerDecoder)

set_target_properties(kmerDecoderTests PROPERTIES CXX_STANDARD 17)

# enable compiler warnings
if(NOT TEST_INSTALLED_VERSION)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(kmerDecoder PUBLIC -Wall -pedantic -Wextra)
    elseif(MSVC)
        target_compile_options(kmerDecoder PUBLIC /W4 /WX)
        target_compile_definitions(kmerDecoderTests PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
    endif()
endif()

enable_testing()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(kmerDecoderTests)

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
    target_compile_options(kmerDecoder PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
    target_link_options(kmerDecoder PUBLIC -fprofile-arcs -ftest-coverage)
endif()