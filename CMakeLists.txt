cmake_minimum_required(VERSION 3.21)

# the run-length encoded BWT, SA samples at the run boundaries and LCP values at the run starts, computed with
# the prefix-free parsing
project(pfp-rlbwt LANGUAGES C CXX)

# add external/string-utils, if the parent project has not added it
if(NOT TARGET string_utils)
  file(GLOB pfp_rlbwt_string_utils_contents "${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/*")

  if(pfp_rlbwt_string_utils_contents STREQUAL "")
    find_package(Git REQUIRED)
    message(STATUS "Downloading the git submodule external/string-utils")
    execute_process(COMMAND "${GIT_EXECUTABLE}" submodule update --init -- external/string-utils
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  # the prefix of the options and environment variables (PFP_RLBWT_MARCH_NATIVE, PFP_RLBWT_TEST_SECONDS, ...)
  if(PROJECT_IS_TOP_LEVEL)
    set(STRING_UTILS_PREFIX PFP_RLBWT)
    include("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/cmake/settings.cmake")
  endif()

  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils" EXCLUDE_FROM_ALL)
endif()

# the libraries of string-utils beyond its core that are used here
string_utils_add_components(libsais)

option(PFP_RLBWT_BUILD_CLI "Build pfp-rlbwt" ${PROJECT_IS_TOP_LEVEL})
option(PFP_RLBWT_BUILD_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
# the compressor (library pfp_rlbwt_comp, CLI pfp-rlbwt-comp and its tests) needs external/move-data-structure, which
# is downloaded only if one of them is built
option(PFP_RLBWT_BUILD_COMP "Build the compressor library pfp_rlbwt_comp" ${PROJECT_IS_TOP_LEVEL})
option(PFP_RLBWT_BUILD_COMP_CLI "Build pfp-rlbwt-comp" ${PROJECT_IS_TOP_LEVEL})

# the libraries the instantiations are compiled and linked with
set(pfp_rlbwt_dependencies libsais libsais40 ips4o ips2ra sux sdsl_core OpenMP::OpenMP_CXX)

# adds the object libraries <prefix>_<bits> with the explicit instantiations (see the end of
# include/pfp_rlbwt/pfp_rlbwt.hpp) for each symbol width and stores their objects in the variable objects; they are
# compiled without the sources string-utils adds to its users (malloc_count), which belong to the programs
function(pfp_rlbwt_add_instantiations prefix objects)
  set(list "")

  foreach(bits 8 16 24 32 40 48)
    add_library(${prefix}_${bits} OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/src/instantiations.cpp")
    target_compile_definitions(${prefix}_${bits} PRIVATE PFP_RLBWT_SYMBOL_BITS=${bits} ${ARGN}
      $<TARGET_PROPERTY:string_utils,INTERFACE_COMPILE_DEFINITIONS>
      $<TARGET_PROPERTY:malloc_count,INTERFACE_COMPILE_DEFINITIONS>)
    target_include_directories(${prefix}_${bits} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/"
      $<TARGET_PROPERTY:string_utils,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:malloc_count,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(${prefix}_${bits} PRIVATE ${pfp_rlbwt_dependencies})
    list(APPEND list $<TARGET_OBJECTS:${prefix}_${bits}>)
  endforeach()

  set(${objects} ${list} PARENT_SCOPE)
endfunction()

pfp_rlbwt_add_instantiations(pfp_rlbwt pfp_rlbwt_objects)
add_library(pfp_rlbwt STATIC ${pfp_rlbwt_objects})
target_include_directories(pfp_rlbwt PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(pfp_rlbwt PUBLIC libsais libsais40 ips4o ips2ra OpenMP::OpenMP_CXX INTERFACE string_utils)

if(PFP_RLBWT_BUILD_CLI)
  string_utils_add_executable(pfp-rlbwt cli/pfp-rlbwt.cpp "${CMAKE_BINARY_DIR}/cli/" pfp_rlbwt)
endif()

if(PFP_RLBWT_BUILD_COMP OR PFP_RLBWT_BUILD_COMP_CLI)
  # add external/move-data-structure, if the parent project has not added it (it uses the string-utils added above)
  if(NOT TARGET move_data_structure)
    string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" PATHS external/move-data-structure)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/move-data-structure" EXCLUDE_FROM_ALL)
  endif()

  add_library(pfp_rlbwt_comp INTERFACE)
  target_link_libraries(pfp_rlbwt_comp INTERFACE pfp_rlbwt move_data_structure string_utils)
endif()

if(PFP_RLBWT_BUILD_COMP_CLI)
  string_utils_add_executable(pfp-rlbwt-comp cli/pfp-rlbwt-comp.cpp "${CMAKE_BINARY_DIR}/cli/" pfp_rlbwt_comp)
endif()

if(PFP_RLBWT_BUILD_TESTS)
  include(GoogleTest)
  enable_testing()
  string_utils_add_components(tests gtl)

  string_utils_add_executable(test-pfp-rlbwt tests/test-pfp-rlbwt.cpp "${CMAKE_BINARY_DIR}/tests/" pfp_rlbwt
    string_utils_test gtl)

  # the same tests, with 40-bit indices for inputs of libsais of more than 2^10 values and 64-bit indices for
  # more than 2^12 values
  set(pfp_rlbwt_wide_definitions PFP_RLBWT_MAX_INDEX32=1024 PFP_RLBWT_MAX_INDEX40=4096)
  pfp_rlbwt_add_instantiations(pfp_rlbwt_wide pfp_rlbwt_wide_objects ${pfp_rlbwt_wide_definitions})
  string_utils_add_executable(test-pfp-rlbwt-wide tests/test-pfp-rlbwt.cpp "${CMAKE_BINARY_DIR}/tests/" string_utils_test
    gtl ${pfp_rlbwt_dependencies})
  target_sources(test-pfp-rlbwt-wide PRIVATE ${pfp_rlbwt_wide_objects})
  target_include_directories(test-pfp-rlbwt-wide PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/")
  target_compile_definitions(test-pfp-rlbwt-wide PRIVATE ${pfp_rlbwt_wide_definitions})

  gtest_discover_tests(test-pfp-rlbwt)
  gtest_discover_tests(test-pfp-rlbwt-wide)

  if(TARGET pfp_rlbwt_comp)
    string_utils_add_executable(test-pfp-rlbwt-comp tests/test-pfp-rlbwt-comp.cpp "${CMAKE_BINARY_DIR}/tests/"
      pfp_rlbwt_comp string_utils_test)
    set_target_properties(test-pfp-rlbwt-comp PROPERTIES INTERPROCEDURAL_OPTIMIZATION OFF)
    gtest_discover_tests(test-pfp-rlbwt-comp)
  endif()
endif()
