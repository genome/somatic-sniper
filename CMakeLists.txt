cmake_minimum_required(VERSION 2.6)

project(somatic_sniper)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)

# TODO: replace with FindSamtools.cmake
set(SAMTOOLS_INC /gsc/pkg/bio/samtools/samtools-0.1.6)
set(SAMTOOLS_LIB /gsc/pkg/bio/samtools/samtools-0.1.6)
include_directories(${SAMTOOLS_INC})
link_directories(${SAMTOOLS_LIB})

# perform out of source build in build directory
add_subdirectory(src build/src)

# unit tests
find_package(GTest)
if(GTEST_FOUND)
    message("Google Test framework found, building unit tests")
    enable_testing(true)
    include_directories(src)
    add_subdirectory(test build/test)
else()
    message("Google Test framework not found, no tests will be built")
    message("GTest is available at http://code.google.com/p/googletest/")
    message("Ubuntu users can likely sudo apt-get install gtest")
endif()
